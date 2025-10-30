"""
qpcr_analysis_app.py — Complete qPCR Analysis Web Application
------------------------------------------------------------
A Streamlit app for qPCR ΔΔCt analysis with:
  - Multiple file upload
  - Automatic data parsing
  - ΔΔCt calculations
  - Statistical analysis
  - Graph generation
  - Results download

Usage:
    streamlit run qpcr_analysis_app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import io
from datetime import datetime
import zipfile

# Set page config
st.set_page_config(
    page_title="qPCR Analysis Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ------------------------------------------------------------
# DATA PARSING FUNCTIONS (Embedded)
# ------------------------------------------------------------

def normalize_header(header: str) -> str:
    """Normalize column headers"""
    import re
    import unicodedata
    if not isinstance(header, str):
        header = str(header)
    header = unicodedata.normalize("NFKD", header).lower().strip()
    header = header.replace("\u0442", "t")
    header = re.sub(r"[^0-9a-z\s]", " ", header)
    return " ".join(header.split())


def find_column(columns, key):
    """Find best matching column name"""
    import difflib
    
    CANONICAL_KEYS = {
        "well": ["well", "position"],
        "sample_name": ["sample name", "sample", "name"],
        "target_name": ["target name", "target", "gene", "reporter", "assay"],
        "ct": ["ct", "cq", "cт", "cq mean", "ct mean"]
    }
    
    candidates = CANONICAL_KEYS[key]
    normalized = {col: normalize_header(col) for col in columns}
    
    for orig, norm in normalized.items():
        if any(cand in norm for cand in candidates):
            return orig
    
    best, score = None, 0
    for orig, norm in normalized.items():
        match = difflib.get_close_matches(norm, candidates, n=1, cutoff=0.0)
        if match:
            s = difflib.SequenceMatcher(a=norm, b=match[0]).ratio()
            if s > score:
                score, best = s, orig
    return best


def to_ct(value):
    """Convert Ct values to numeric"""
    import re
    if pd.isna(value):
        return np.nan
    s = str(value).strip()
    if s.lower() in ["undetermined", "undet", "na", "n/a", "-", ""]:
        return np.nan
    s = re.sub(r"[^0-9eE\.\-]", "", s)
    try:
        return float(s)
    except ValueError:
        return np.nan


def parse_qpcr_file(uploaded_file):
    """Parse uploaded qPCR file"""
    try:
        file_name = uploaded_file.name
        
        if file_name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
            sheets = {"csv": df}
        else:
            sheets = pd.read_excel(uploaded_file, sheet_name=None, engine="openpyxl")
        
        parsed_list = []
        
        for sheet_name, df in sheets.items():
            if df.empty:
                continue
            
            # Check if first row is header
            if not all(isinstance(c, str) for c in df.columns):
                potential_header = df.iloc[0]
                if potential_header.isna().sum() < len(potential_header) - 2:
                    df.columns = potential_header
                    df = df.drop(df.index[0]).reset_index(drop=True)
            
            cols = [str(c) for c in df.columns]
            col_well = find_column(cols, "well")
            col_sample = find_column(cols, "sample_name")
            col_target = find_column(cols, "target_name")
            col_ct = find_column(cols, "ct")
            
            if not col_ct:
                continue
            
            parsed = pd.DataFrame()
            parsed["ct"] = df[col_ct].apply(to_ct)
            parsed["well"] = df[col_well] if col_well else np.nan
            parsed["sample_name"] = df[col_sample].astype(str).str.strip() if col_sample else parsed["well"]
            parsed["target_name"] = df[col_target].astype(str).str.strip().str.upper() if col_target else np.nan
            parsed["file_source"] = file_name
            parsed["sheet_name"] = sheet_name
            
            parsed = parsed.dropna(subset=["ct"], how="all")
            parsed = parsed[["well", "sample_name", "target_name", "ct", "file_source", "sheet_name"]]
            
            if not parsed.empty:
                parsed_list.append(parsed)
        
        if parsed_list:
            return pd.concat(parsed_list, ignore_index=True)
        else:
            return pd.DataFrame()
            
    except Exception as e:
        st.error(f"Error parsing {uploaded_file.name}: {str(e)}")
        return pd.DataFrame()


# ------------------------------------------------------------
# CALCULATION FUNCTIONS (Embedded)
# ------------------------------------------------------------

class QPCRCalculator:
    """qPCR ΔΔCt Calculator"""
    
    def __init__(self, data, reference_gene, control_group):
        self.data = data.copy()
        self.reference_gene = reference_gene.upper()
        self.control_group = control_group
        self.validate_data()
    
    def validate_data(self):
        """Validate input data"""
        required_cols = ["sample_name", "target_name", "ct"]
        missing = [col for col in required_cols if col not in self.data.columns]
        
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        if self.reference_gene not in self.data["target_name"].unique():
            available = self.data["target_name"].unique()
            raise ValueError(
                f"Reference gene '{self.reference_gene}' not found.\n"
                f"Available genes: {list(available)}"
            )
    
    def calculate_all(self):
        """Run complete analysis"""
        # Step 1: Replicate statistics
        grouped = self.data.groupby(["sample_name", "target_name"])["ct"].agg([
            ("ct_mean", "mean"),
            ("ct_sd", "std"),
            ("ct_count", "count")
        ]).reset_index()
        
        grouped["ct_cv_percent"] = (grouped["ct_sd"] / grouped["ct_mean"]) * 100
        grouped["high_variation"] = grouped["ct_cv_percent"] > 5.0
        self.replicate_stats = grouped
        
        # Step 2: ΔCt calculation
        ref_data = grouped[grouped["target_name"] == self.reference_gene][
            ["sample_name", "ct_mean"]
        ].rename(columns={"ct_mean": "ref_ct"})
        
        target_data = grouped[grouped["target_name"] != self.reference_gene].copy()
        merged = target_data.merge(ref_data, on="sample_name", how="left")
        merged["delta_ct"] = merged["ct_mean"] - merged["ref_ct"]
        
        ref_sd = grouped[grouped["target_name"] == self.reference_gene][
            ["sample_name", "ct_sd"]
        ].rename(columns={"ct_sd": "ref_sd"})
        
        merged = merged.merge(ref_sd, on="sample_name", how="left")
        merged["delta_ct_sd"] = np.sqrt(merged["ct_sd"]**2 + merged["ref_sd"]**2)
        self.delta_ct_data = merged
        
        # Step 3: ΔΔCt and Fold Change
        control_data = merged[
            merged["sample_name"].str.contains(self.control_group, case=False)
        ].groupby("target_name")["delta_ct"].mean().reset_index()
        control_data.columns = ["target_name", "control_delta_ct"]
        
        merged = merged.merge(control_data, on="target_name", how="left")
        merged["delta_delta_ct"] = merged["delta_ct"] - merged["control_delta_ct"]
        merged["fold_change"] = 2 ** (-merged["delta_delta_ct"])
        merged["fold_change_sd"] = merged["fold_change"] * np.log(2) * merged["delta_ct_sd"]
        self.final_results = merged
        
        # Step 4: Statistics
        stats_results = []
        for target in merged["target_name"].unique():
            target_data = merged[merged["target_name"] == target]
            
            control_values = target_data[
                target_data["sample_name"].str.contains(self.control_group, case=False)
            ]["delta_ct"].values
            
            for group in target_data["sample_name"].unique():
                if self.control_group.lower() in group.lower():
                    continue
                
                treatment_values = target_data[
                    target_data["sample_name"] == group
                ]["delta_ct"].values
                
                if len(treatment_values) >= 2 and len(control_values) >= 2:
                    t_stat, p_value = stats.ttest_ind(treatment_values, control_values)
                else:
                    p_value = np.nan
                
                if pd.isna(p_value):
                    significance = "N/A"
                elif p_value < 0.001:
                    significance = "***"
                elif p_value < 0.01:
                    significance = "**"
                elif p_value < 0.05:
                    significance = "*"
                else:
                    significance = "ns"
                
                stats_results.append({
                    "target_name": target,
                    "sample_name": group,
                    "p_value": p_value,
                    "significance": significance
                })
        
        self.statistics = pd.DataFrame(stats_results)
        self.final_results = self.final_results.merge(
            self.statistics,
            on=["target_name", "sample_name"],
            how="left"
        )
        
        return self.get_summary()
    
    def get_summary(self):
        """Get summary table"""
        summary = self.final_results[[
            "sample_name", "target_name", "ct_mean", "ct_sd",
            "delta_ct", "fold_change", "fold_change_sd", "p_value", "significance"
        ]].copy()
        
        numeric_cols = ["ct_mean", "ct_sd", "delta_ct", "fold_change", "fold_change_sd", "p_value"]
        summary[numeric_cols] = summary[numeric_cols].round(3)
        
        return summary


# ------------------------------------------------------------
# VISUALIZATION FUNCTIONS (Embedded)
# ------------------------------------------------------------

def create_bar_graph(data, gene, style="default"):
    """Create bar chart for one gene"""
    gene_data = data[data["target_name"].str.upper() == gene.upper()].copy()
    
    if gene_data.empty:
        return None
    
    gene_data = gene_data.sort_values("sample_name")
    
    samples = gene_data["sample_name"].values
    fold_changes = gene_data["fold_change"].values
    errors = gene_data["fold_change_sd"].values
    sig_markers = gene_data["significance"].values
    
    # Color scheme
    colors = []
    for sample in samples:
        if any(kw in sample.lower() for kw in ["control", "ctrl", "vehicle", "dmso"]):
            colors.append("#95a5a6")
        else:
            colors.append("#3498db")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    x_pos = np.arange(len(samples))
    bars = ax.bar(x_pos, fold_changes, color=colors, alpha=0.8, edgecolor="black", linewidth=1.2)
    
    ax.errorbar(x_pos, fold_changes, yerr=errors, fmt='none', ecolor='black', 
                capsize=5, capthick=1.5, linewidth=1.5)
    
    # Add significance markers
    for i, (fc, err, sig) in enumerate(zip(fold_changes, errors, sig_markers)):
        if pd.notna(sig) and sig in ['*', '**', '***']:
            y_pos = fc + err + max(fold_changes) * 0.05
            ax.text(i, y_pos, sig, ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.set_ylabel("Relative mRNA Expression (Fold Change)", fontsize=12, fontweight='bold')
    ax.set_xlabel("Sample", fontsize=12, fontweight='bold')
    ax.set_title(f"{gene} Expression", fontsize=14, fontweight='bold', pad=20)
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_ylim(bottom=0, top=max(fold_changes + errors) * 1.2)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    return fig


# ------------------------------------------------------------
# STREAMLIT APP
# ------------------------------------------------------------

def main():
    
    # Title and description
    st.title("🧬 qPCR ΔΔCt Analysis Tool")
    st.markdown("""
    Upload your qPCR raw data files and get automated analysis with:
    - ΔΔCt calculations and fold change
    - Statistical analysis (t-tests)
    - Publication-ready graphs
    - Downloadable Excel results
    """)
    
    # Sidebar configuration
    st.sidebar.header("⚙️ Configuration")
    
    reference_gene = st.sidebar.text_input(
        "Reference Gene (Housekeeping)",
        value="ACTB",
        help="Enter the gene symbol for your internal control (e.g., ACTB, GAPDH, 18S)"
    )
    
    control_group = st.sidebar.text_input(
        "Control Group Identifier",
        value="Control",
        help="Text in sample names that identifies control samples (e.g., Control, Vehicle, DMSO)"
    )
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    ### 📋 Instructions
    1. Upload one or more Excel/CSV files
    2. Configure reference gene and control group
    3. Click 'Analyze Data'
    4. Download results and graphs
    
    ### 📝 File Format
    Your files should contain columns:
    - Well (optional)
    - Sample Name
    - Target Name (Gene)
    - Ct value
    """)
    
    # Main content
    st.header("📤 Step 1: Upload Data Files")
    
    uploaded_files = st.file_uploader(
        "Upload qPCR Excel or CSV files",
        type=["xlsx", "xls", "csv"],
        accept_multiple_files=True,
        help="You can upload multiple files at once"
    )
    
    if uploaded_files:
        st.success(f"✅ Uploaded {len(uploaded_files)} file(s)")
        
        # Show uploaded files
        with st.expander("📁 Uploaded Files"):
            for file in uploaded_files:
                st.write(f"• {file.name}")
        
        # Parse files button
        if st.button("🔍 Parse and Preview Data", type="primary"):
            with st.spinner("Parsing files..."):
                
                # Parse all files
                all_parsed = []
                parse_errors = []
                
                for file in uploaded_files:
                    parsed = parse_qpcr_file(file)
                    if not parsed.empty:
                        all_parsed.append(parsed)
                    else:
                        parse_errors.append(file.name)
                
                if parse_errors:
                    st.warning(f"⚠️ Could not parse: {', '.join(parse_errors)}")
                
                if all_parsed:
                    combined_data = pd.concat(all_parsed, ignore_index=True)
                    st.session_state['parsed_data'] = combined_data
                    
                    st.success(f"✅ Successfully parsed {len(combined_data)} rows")
                    
                    # Show preview
                    st.subheader("📊 Data Preview")
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Rows", len(combined_data))
                    with col2:
                        st.metric("Unique Samples", combined_data['sample_name'].nunique())
                    with col3:
                        st.metric("Unique Genes", combined_data['target_name'].nunique())
                    
                    # Show data table
                    st.dataframe(combined_data.head(20), use_container_width=True)
                    
                    # Show genes found
                    st.markdown("**🧬 Genes Found:**")
                    genes = combined_data['target_name'].unique()
                    st.write(", ".join(sorted(genes)))
                    
                    # Show samples found
                    st.markdown("**🔬 Samples Found:**")
                    samples = combined_data['sample_name'].unique()
                    st.write(", ".join(sorted(samples)))
                    
                else:
                    st.error("❌ No valid data found in uploaded files")
    
    # Analysis section
    if 'parsed_data' in st.session_state:
        st.markdown("---")
        st.header("🧮 Step 2: Run Analysis")
        
        parsed_data = st.session_state['parsed_data']
        
        # Validate reference gene
        available_genes = parsed_data['target_name'].unique()
        if reference_gene.upper() not in available_genes:
            st.error(f"⚠️ Reference gene '{reference_gene}' not found in data!")
            st.info(f"Available genes: {', '.join(available_genes)}")
            st.stop()
        
        if st.button("🚀 Analyze Data", type="primary"):
            with st.spinner("Running ΔΔCt analysis..."):
                
                try:
                    # Run calculation
                    calculator = QPCRCalculator(
                        data=parsed_data,
                        reference_gene=reference_gene,
                        control_group=control_group
                    )
                    
                    summary = calculator.calculate_all()
                    
                    st.session_state['calculator'] = calculator
                    st.session_state['summary'] = summary
                    
                    st.success("✅ Analysis complete!")
                    
                    # Display results
                    st.subheader("📊 Results Summary")
                    
                    # Metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        sig_count = summary[summary['significance'].isin(['*', '**', '***'])].shape[0]
                        st.metric("Significant Results", sig_count)
                    with col2:
                        avg_fc = summary['fold_change'].mean()
                        st.metric("Avg Fold Change", f"{avg_fc:.2f}")
                    with col3:
                        high_cv = calculator.replicate_stats[
                            calculator.replicate_stats['ct_cv_percent'] > 5
                        ].shape[0]
                        st.metric("High CV% Samples", high_cv)
                    
                    # Show summary table
                    st.dataframe(summary, use_container_width=True)
                    
                    # Show significance breakdown
                    st.markdown("**📈 Statistical Significance:**")
                    sig_counts = summary['significance'].value_counts()
                    for sig, count in sig_counts.items():
                        st.write(f"• {sig}: {count}")
                    
                except Exception as e:
                    st.error(f"❌ Analysis failed: {str(e)}")
                    import traceback
                    st.code(traceback.format_exc())
    
    # Download section
    if 'calculator' in st.session_state and 'summary' in st.session_state:
        st.markdown("---")
        st.header("📥 Step 3: Download Results")
        
        calculator = st.session_state['calculator']
        summary = st.session_state['summary']
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("📊 Excel Results")
            
            # Create Excel file in memory
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                calculator.replicate_stats.to_excel(writer, sheet_name='Replicate_Stats', index=False)
                calculator.delta_ct_data.to_excel(writer, sheet_name='Delta_Ct', index=False)
                calculator.final_results.to_excel(writer, sheet_name='Final_Results', index=False)
                summary.to_excel(writer, sheet_name='Summary', index=False)
            
            excel_data = output.getvalue()
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            st.download_button(
                label="📥 Download Excel Results",
                data=excel_data,
                file_name=f"qpcr_results_{timestamp}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        with col2:
            st.subheader("📈 Generate Graphs")
            
            genes = summary['target_name'].unique()
            selected_genes = st.multiselect(
                "Select genes to plot",
                options=genes,
                default=list(genes)[:3] if len(genes) >= 3 else list(genes)
            )
            
            if st.button("Generate Graphs"):
                if selected_genes:
                    with st.spinner("Creating graphs..."):
                        
                        # Create ZIP file with all graphs
                        zip_buffer = io.BytesIO()
                        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                            
                            for gene in selected_genes:
                                fig = create_bar_graph(summary, gene)
                                if fig:
                                    # Save figure to buffer
                                    img_buffer = io.BytesIO()
                                    fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                                    img_buffer.seek(0)
                                    
                                    # Add to ZIP
                                    zip_file.writestr(f"{gene}_expression.png", img_buffer.getvalue())
                                    
                                    # Display in Streamlit
                                    st.pyplot(fig)
                                    plt.close(fig)
                        
                        zip_buffer.seek(0)
                        
                        st.download_button(
                            label="📥 Download All Graphs (ZIP)",
                            data=zip_buffer.getvalue(),
                            file_name=f"qpcr_graphs_{timestamp}.zip",
                            mime="application/zip"
                        )
                else:
                    st.warning("Please select at least one gene to plot")


if __name__ == "__main__":
    main()