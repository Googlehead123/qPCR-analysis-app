"""
qpcr_analysis_app.py — Complete qPCR Analysis Web Application
------------------------------------------------------------
A Streamlit app for qPCR ΔΔCt analysis with:
  - Multiple file upload
  - Parse only: Well, Sample Name, Target Name, Ct
  - User-configurable reference group and control gene
  - Relative gene expression calculations
  - Statistical analysis
  - Excel and graph downloads

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
import re
import unicodedata

# Set page config
st.set_page_config(
    page_title="qPCR Analysis Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ------------------------------------------------------------
# DATA PARSING FUNCTIONS
# ------------------------------------------------------------

def normalize_header(header: str) -> str:
    """Normalize column headers to lowercase and remove special chars"""
    if not isinstance(header, str):
        header = str(header)
    # Normalize unicode
    header = unicodedata.normalize("NFKD", header).lower().strip()
    # Replace Cyrillic т with Latin t
    header = header.replace("\u0442", "t")
    # Remove special characters except spaces
    header = re.sub(r"[^0-9a-z\s]", " ", header)
    return " ".join(header.split())


def find_ct_column(columns):
    """Find Ct column (handles Ct, CT, Cт, Cq variations)"""
    ct_patterns = ["ct", "cт", "cq"]
    
    for col in columns:
        normalized = normalize_header(col)
        # Check for exact match or contains pattern
        if normalized in ct_patterns or any(pattern in normalized for pattern in ct_patterns):
            return col
    return None


def find_column_fuzzy(columns, target_words):
    """Find column by matching target words"""
    for col in columns:
        normalized = normalize_header(col)
        if any(word in normalized for word in target_words):
            return col
    return None


def parse_qpcr_file(uploaded_file):
    """
    Parse uploaded qPCR file - Extract ONLY 4 columns:
    Well, Sample Name, Target Name, Ct
    """
    try:
        file_name = uploaded_file.name
        
        # Read file
        if file_name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
            sheets = {"csv": df}
        else:
            sheets = pd.read_excel(uploaded_file, sheet_name=None, engine="openpyxl")
        
        parsed_list = []
        
        for sheet_name, df in sheets.items():
            if df.empty:
                continue
            
            # Check if first row should be header
            if not all(isinstance(c, str) for c in df.columns):
                potential_header = df.iloc[0]
                if potential_header.notna().sum() >= 3:  # At least 3 non-null values
                    df.columns = potential_header
                    df = df.drop(df.index[0]).reset_index(drop=True)
            
            # Find required columns
            cols = [str(c) for c in df.columns]
            
            # 1. Find Well column
            col_well = find_column_fuzzy(cols, ["well", "position"])
            
            # 2. Find Sample Name column
            col_sample = find_column_fuzzy(cols, ["sample name", "sample", "name"])
            
            # 3. Find Target Name column
            col_target = find_column_fuzzy(cols, ["target name", "target", "gene", "reporter", "assay"])
            
            # 4. Find Ct column (CRITICAL - handles Ct, CT, Cт, Cq)
            col_ct = find_ct_column(cols)
            
            # Skip if essential columns not found
            if not col_ct:
                st.warning(f"⚠️ No Ct column found in {file_name} - {sheet_name}")
                continue
            
            if not col_sample:
                st.warning(f"⚠️ No Sample Name column found in {file_name} - {sheet_name}")
                continue
                
            if not col_target:
                st.warning(f"⚠️ No Target Name column found in {file_name} - {sheet_name}")
                continue
            
            # Extract only the 4 required columns
            parsed = pd.DataFrame()
            
            # Well (optional - if not found, use row number)
            if col_well:
                parsed["Well"] = df[col_well]
            else:
                parsed["Well"] = [f"Row_{i+1}" for i in range(len(df))]
            
            # Sample Name
            parsed["Sample Name"] = df[col_sample].astype(str).str.strip()
            
            # Target Name
            parsed["Target Name"] = df[col_target].astype(str).str.strip().str.upper()
            
            # Ct - Convert to numeric, handle "Undetermined"
            def parse_ct(value):
                if pd.isna(value):
                    return np.nan
                s = str(value).strip().lower()
                if s in ["undetermined", "undet", "na", "n/a", "-", ""]:
                    return np.nan
                # Remove non-numeric characters except decimal point and minus
                s = re.sub(r"[^0-9\.\-]", "", s)
                try:
                    return float(s)
                except ValueError:
                    return np.nan
            
            parsed["Ct"] = df[col_ct].apply(parse_ct)
            
            # Add metadata
            parsed["File Source"] = file_name
            parsed["Sheet Name"] = sheet_name
            
            # Remove rows where Ct is NaN
            parsed = parsed.dropna(subset=["Ct"])
            
            # Remove empty sample names
            parsed = parsed[parsed["Sample Name"].str.len() > 0]
            parsed = parsed[parsed["Target Name"].str.len() > 0]
            
            if not parsed.empty:
                parsed_list.append(parsed)
        
        if parsed_list:
            combined = pd.concat(parsed_list, ignore_index=True)
            return combined
        else:
            return pd.DataFrame()
            
    except Exception as e:
        st.error(f"Error parsing {uploaded_file.name}: {str(e)}")
        return pd.DataFrame()


# ------------------------------------------------------------
# CALCULATION FUNCTIONS
# ------------------------------------------------------------

class QPCRCalculator:
    """qPCR ΔΔCt Calculator with user-defined reference group and control gene"""
    
    def __init__(self, data, control_gene, reference_group):
        """
        Initialize calculator
        
        Args:
            data: DataFrame with columns [Well, Sample Name, Target Name, Ct]
            control_gene: Reference/housekeeping gene (e.g., ACTB, GAPDH)
            reference_group: Control/reference group in Sample Name (e.g., Control, Vehicle)
        """
        self.data = data.copy()
        self.control_gene = control_gene.upper()
        self.reference_group = reference_group
        self.validate_data()
    
    def validate_data(self):
        """Validate input data"""
        required_cols = ["Sample Name", "Target Name", "Ct"]
        missing = [col for col in required_cols if col not in self.data.columns]
        
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        # Check if control gene exists
        available_genes = self.data["Target Name"].unique()
        if self.control_gene not in available_genes:
            raise ValueError(
                f"Control gene '{self.control_gene}' not found in data.\n"
                f"Available genes: {list(available_genes)}"
            )
        
        # Check if reference group exists
        available_samples = self.data["Sample Name"].unique()
        matching_samples = [s for s in available_samples if self.reference_group.lower() in s.lower()]
        
        if not matching_samples:
            raise ValueError(
                f"Reference group '{self.reference_group}' not found in Sample Names.\n"
                f"Available samples: {list(available_samples)}"
            )
    
    def calculate_all(self):
        """Run complete ΔΔCt analysis"""
        
        # STEP 1: Calculate mean Ct for technical replicates
        st.write("### 📊 Step 1: Calculating replicate statistics...")
        
        grouped = self.data.groupby(["Sample Name", "Target Name"])["Ct"].agg([
            ("Ct_Mean", "mean"),
            ("Ct_SD", "std"),
            ("Ct_Count", "count")
        ]).reset_index()
        
        # Calculate CV% (Coefficient of Variation)
        grouped["CV_Percent"] = (grouped["Ct_SD"] / grouped["Ct_Mean"]) * 100
        grouped["High_Variation"] = grouped["CV_Percent"] > 5.0
        
        self.replicate_stats = grouped
        
        # Warn about high variation
        high_var = grouped[grouped["High_Variation"]]
        if not high_var.empty:
            st.warning(f"⚠️ Found {len(high_var)} sample(s) with CV > 5%")
            with st.expander("View samples with high variation"):
                st.dataframe(high_var[["Sample Name", "Target Name", "Ct_Mean", "CV_Percent"]])
        
        # STEP 2: Calculate ΔCt (Target - Control Gene)
        st.write("### 📐 Step 2: Calculating ΔCt values...")
        
        # Separate control gene data
        control_gene_data = grouped[grouped["Target Name"] == self.control_gene][
            ["Sample Name", "Ct_Mean"]
        ].rename(columns={"Ct_Mean": "Control_Gene_Ct"})
        
        # Merge with target genes
        target_data = grouped[grouped["Target Name"] != self.control_gene].copy()
        merged = target_data.merge(control_gene_data, on="Sample Name", how="left")
        
        # Calculate ΔCt = Ct(target) - Ct(control gene)
        merged["Delta_Ct"] = merged["Ct_Mean"] - merged["Control_Gene_Ct"]
        
        # Error propagation: SD(ΔCt) = sqrt(SD(target)² + SD(control)²)
        control_sd = grouped[grouped["Target Name"] == self.control_gene][
            ["Sample Name", "Ct_SD"]
        ].rename(columns={"Ct_SD": "Control_Gene_SD"})
        
        merged = merged.merge(control_sd, on="Sample Name", how="left")
        merged["Delta_Ct_SD"] = np.sqrt(merged["Ct_SD"]**2 + merged["Control_Gene_SD"]**2)
        
        self.delta_ct_data = merged
        
        # STEP 3: Calculate ΔΔCt (Treatment - Reference Group)
        st.write("### 🔬 Step 3: Calculating ΔΔCt and Relative Expression...")
        
        # Get reference group ΔCt values
        reference_data = merged[
            merged["Sample Name"].str.contains(self.reference_group, case=False)
        ].groupby("Target Name")["Delta_Ct"].mean().reset_index()
        reference_data.columns = ["Target Name", "Reference_Delta_Ct"]
        
        # Merge with all data
        merged = merged.merge(reference_data, on="Target Name", how="left")
        
        # Calculate ΔΔCt = ΔCt(sample) - ΔCt(reference)
        merged["Delta_Delta_Ct"] = merged["Delta_Ct"] - merged["Reference_Delta_Ct"]
        
        # Calculate Relative Expression (Fold Change) = 2^(-ΔΔCt)
        merged["Relative_Expression"] = 2 ** (-merged["Delta_Delta_Ct"])
        
        # Error propagation for relative expression
        merged["Relative_Expression_SD"] = merged["Relative_Expression"] * np.log(2) * merged["Delta_Ct_SD"]
        
        self.final_results = merged
        
        # STEP 4: Statistical Analysis (t-test)
        st.write("### 📈 Step 4: Performing statistical analysis...")
        
        stats_results = []
        
        for target in merged["Target Name"].unique():
            target_data = merged[merged["Target Name"] == target]
            
            # Get reference group Delta_Ct values
            reference_values = target_data[
                target_data["Sample Name"].str.contains(self.reference_group, case=False)
            ]["Delta_Ct"].values
            
            # Compare each non-reference group to reference
            for sample in target_data["Sample Name"].unique():
                if self.reference_group.lower() in sample.lower():
                    # This is the reference group - set p=1, significance=ns
                    stats_results.append({
                        "Target Name": target,
                        "Sample Name": sample,
                        "P_Value": 1.0,
                        "Significance": "Reference"
                    })
                    continue
                
                # Get treatment values
                treatment_values = target_data[
                    target_data["Sample Name"] == sample
                ]["Delta_Ct"].values
                
                # Perform t-test if enough replicates
                if len(treatment_values) >= 2 and len(reference_values) >= 2:
                    try:
                        t_stat, p_value = stats.ttest_ind(treatment_values, reference_values)
                    except:
                        p_value = np.nan
                else:
                    p_value = np.nan
                
                # Determine significance level
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
                    "Target Name": target,
                    "Sample Name": sample,
                    "P_Value": p_value,
                    "Significance": significance
                })
        
        self.statistics = pd.DataFrame(stats_results)
        
        # Merge statistics with final results
        self.final_results = self.final_results.merge(
            self.statistics,
            on=["Target Name", "Sample Name"],
            how="left"
        )
        
        return self.get_summary()
    
    def get_summary(self):
        """Generate clean summary table"""
        summary = self.final_results[[
            "Sample Name",
            "Target Name",
            "Ct_Mean",
            "Ct_SD",
            "Delta_Ct",
            "Relative_Expression",
            "Relative_Expression_SD",
            "P_Value",
            "Significance"
        ]].copy()
        
        # Round numeric columns
        numeric_cols = ["Ct_Mean", "Ct_SD", "Delta_Ct", "Relative_Expression", "Relative_Expression_SD", "P_Value"]
        for col in numeric_cols:
            if col in summary.columns:
                summary[col] = summary[col].round(3)
        
        return summary


# ------------------------------------------------------------
# VISUALIZATION FUNCTIONS
# ------------------------------------------------------------

def create_expression_graph(data, gene):
    """Create bar chart for relative gene expression"""
    gene_data = data[data["Target Name"].str.upper() == gene.upper()].copy()
    
    if gene_data.empty:
        return None
    
    # Sort by sample name
    gene_data = gene_data.sort_values("Sample Name")
    
    samples = gene_data["Sample Name"].values
    expression = gene_data["Relative_Expression"].values
    errors = gene_data["Relative_Expression_SD"].values
    sig_markers = gene_data["Significance"].values
    
    # Assign colors (reference group = gray, others = blue)
    colors = []
    for i, sample in enumerate(samples):
        if sig_markers[i] == "Reference":
            colors.append("#95a5a6")  # Gray for reference
        else:
            colors.append("#3498db")  # Blue for treatments
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x_pos = np.arange(len(samples))
    bars = ax.bar(x_pos, expression, color=colors, alpha=0.8, edgecolor="black", linewidth=1.2)
    
    # Add error bars
    ax.errorbar(x_pos, expression, yerr=errors, fmt='none', ecolor='black',
                capsize=5, capthick=1.5, linewidth=1.5)
    
    # Add significance markers
    for i, (exp, err, sig) in enumerate(zip(expression, errors, sig_markers)):
        if sig in ['*', '**', '***']:
            y_pos = exp + err + max(expression) * 0.05
            ax.text(i, y_pos, sig, ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    # Customize plot
    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel("Relative Gene Expression (Fold Change)", fontsize=12, fontweight='bold')
    ax.set_xlabel("Sample", fontsize=12, fontweight='bold')
    ax.set_title(f"{gene} Relative Expression", fontsize=14, fontweight='bold', pad=20)
    
    # Add reference line at y=1
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.7, label='Reference')
    
    # Set y-axis limits
    y_max = max(expression + errors) * 1.2
    ax.set_ylim(bottom=0, top=y_max)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3)
    
    # Add legend
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    return fig


# ------------------------------------------------------------
# STREAMLIT APP
# ------------------------------------------------------------

def main():
    
    # Title
    st.title("🧬 qPCR Relative Gene Expression Analysis")
    st.markdown("""
    Upload your qPCR data and get automated ΔΔCt analysis with:
    - Relative gene expression calculations (2^-ΔΔCt method)
    - Statistical significance testing
    - Publication-ready graphs
    - Downloadable Excel results
    """)
    
    # Sidebar instructions
    st.sidebar.header("📋 Instructions")
    st.sidebar.markdown("""
    ### How to Use:
    1. **Upload** Excel/CSV files with qPCR data
    2. **Parse** data to extract required columns
    3. **Configure** control gene and reference group
    4. **Analyze** to calculate relative expression
    5. **Download** results and graphs
    
    ### Required Columns:
    - **Well** (optional)
    - **Sample Name** (required)
    - **Target Name** (required)
    - **Ct** or **CT** or **Cт** (required)
    """)
    
    # ========================================
    # STEP 1: Upload Files
    # ========================================
    st.header("📤 Step 1: Upload Data Files")
    
    uploaded_files = st.file_uploader(
        "Upload qPCR Excel or CSV files",
        type=["xlsx", "xls", "csv"],
        accept_multiple_files=True,
        help="Upload one or more files containing qPCR raw data"
    )
    
    if uploaded_files:
        st.success(f"✅ Uploaded {len(uploaded_files)} file(s)")
        
        with st.expander("📁 View uploaded files"):
            for file in uploaded_files:
                st.write(f"• {file.name}")
        
        # Parse button
        if st.button("🔍 Parse Data", type="primary"):
            with st.spinner("Parsing files..."):
                
                all_parsed = []
                
                for file in uploaded_files:
                    parsed = parse_qpcr_file(file)
                    if not parsed.empty:
                        all_parsed.append(parsed)
                
                if all_parsed:
                    combined_data = pd.concat(all_parsed, ignore_index=True)
                    st.session_state['parsed_data'] = combined_data
                    
                    st.success(f"✅ Successfully parsed {len(combined_data)} rows")
                    
                    # Show preview
                    st.subheader("📊 Parsed Data Preview")
                    
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("Total Rows", len(combined_data))
                    with col2:
                        st.metric("Samples", combined_data['Sample Name'].nunique())
                    with col3:
                        st.metric("Genes", combined_data['Target Name'].nunique())
                    with col4:
                        st.metric("Files", combined_data['File Source'].nunique())
                    
                    # Show data table
                    st.dataframe(combined_data[["Well", "Sample Name", "Target Name", "Ct"]].head(20), 
                               use_container_width=True)
                    
                    # Show available genes
                    st.markdown("**🧬 Available Genes:**")
                    genes = sorted(combined_data['Target Name'].unique())
                    st.write(", ".join(genes))
                    
                    # Show available samples
                    st.markdown("**🔬 Available Samples:**")
                    samples = sorted(combined_data['Sample Name'].unique())
                    st.write(", ".join(samples))
                    
                else:
                    st.error("❌ No valid data found in uploaded files")
    
    # ========================================
    # STEP 2: Configuration
    # ========================================
    if 'parsed_data' in st.session_state:
        st.markdown("---")
        st.header("⚙️ Step 2: Configure Analysis Parameters")
        
        parsed_data = st.session_state['parsed_data']
        available_genes = sorted(parsed_data['Target Name'].unique())
        available_samples = sorted(parsed_data['Sample Name'].unique())
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("🧬 Control Gene (Housekeeping)")
            st.markdown("Select the reference/housekeeping gene for normalization")
            
            control_gene = st.selectbox(
                "Control Gene",
                options=available_genes,
                help="Common housekeeping genes: ACTB, GAPDH, 18S, B2M"
            )
            
            st.info(f"✅ Selected: **{control_gene}**")
        
        with col2:
            st.subheader("📊 Reference Group (Control)")
            st.markdown("Enter the identifier for your reference/control group")
            
            # Show sample names as reference
            with st.expander("View all sample names"):
                for sample in available_samples:
                    st.write(f"• {sample}")
            
            reference_group = st.text_input(
                "Reference Group",
                value="Control",
                help="Enter text that identifies control samples (e.g., 'Control', 'Vehicle', 'DMSO')"
            )
            
            # Validate reference group
            matching_samples = [s for s in available_samples if reference_group.lower() in s.lower()]
            if matching_samples:
                st.success(f"✅ Found {len(matching_samples)} matching sample(s):")
                for sample in matching_samples:
                    st.write(f"  • {sample}")
            else:
                st.error("❌ No samples found matching this reference group")
        
        # Save configuration
        if control_gene and reference_group:
            st.session_state['control_gene'] = control_gene
            st.session_state['reference_group'] = reference_group
    
    # ========================================
    # STEP 3: Analysis
    # ========================================
    if 'parsed_data' in st.session_state and 'control_gene' in st.session_state and 'reference_group' in st.session_state:
        st.markdown("---")
        st.header("🚀 Step 3: Run Analysis")
        
        if st.button("🧮 Calculate Relative Expression", type="primary"):
            
            parsed_data = st.session_state['parsed_data']
            control_gene = st.session_state['control_gene']
            reference_group = st.session_state['reference_group']
            
            try:
                with st.spinner("Running ΔΔCt analysis..."):
                    
                    # Initialize calculator
                    calculator = QPCRCalculator(
                        data=parsed_data,
                        control_gene=control_gene,
                        reference_group=reference_group
                    )
                    
                    # Run calculations
                    summary = calculator.calculate_all()
                    
                    # Save to session state
                    st.session_state['calculator'] = calculator
                    st.session_state['summary'] = summary
                    
                    st.success("✅ Analysis Complete!")
                    
                    # Display summary
                    st.subheader("📊 Results Summary")
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        sig_count = summary[summary['Significance'].isin(['*', '**', '***'])].shape[0]
                        st.metric("Significant Results", sig_count)
                    with col2:
                        genes_analyzed = summary['Target Name'].nunique()
                        st.metric("Genes Analyzed", genes_analyzed)
                    with col3:
                        samples_analyzed = summary['Sample Name'].nunique()
                        st.metric("Samples Analyzed", samples_analyzed)
                    
                    # Show summary table
                    st.dataframe(summary, use_container_width=True)
                    
                    # Significance breakdown
                    st.markdown("**📈 Statistical Significance:**")
                    sig_counts = summary['Significance'].value_counts()
                    for sig, count in sig_counts.items():
                        st.write(f"  • {sig}: {count}")
                
            except Exception as e:
                st.error(f"❌ Analysis failed: {str(e)}")
                st.code(str(e))
    
    # ========================================
    # STEP 4: Download Results
    # ========================================
    if 'calculator' in st.session_state and 'summary' in st.session_state:
        st.markdown("---")
        st.header("📥 Step 4: Download Results")
        
        calculator = st.session_state['calculator']
        summary = st.session_state['summary']
        
        col1, col2 = st.columns(2)
        
        # Excel Download
        with col1:
            st.subheader("📊 Excel Results")
            
            # Create Excel file
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
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                help="Download complete analysis results in Excel format"
            )
        
        # Graph Generation
        with col2:
            st.subheader("📈 Generate Graphs")
            
            genes = summary['Target Name'].unique()
            selected_genes = st.multiselect(
                "Select genes to plot",
                options=genes,
                default=list(genes)
            )
            
            if st.button("Generate Graphs"):
                if selected_genes:
                    with st.spinner("Creating graphs..."):
                        
                        # Create ZIP file
                        zip_buffer = io.BytesIO()
                        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                            
                            for gene in selected_genes:
                                fig = create_expression_graph(summary, gene)
                                if fig:
                                    # Save to buffer
                                    img_buffer = io.BytesIO()
                                    fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                                    img_buffer.seek(0)
                                    
                                    # Add to ZIP
                                    zip_file.writestr(f"{gene}_expression.png", img_buffer.getvalue())
                                    
                                    # Display
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