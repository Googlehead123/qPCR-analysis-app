"""
qpcr_analysis_app.py — Complete qPCR Analysis Web Application
------------------------------------------------------------
Features:
  - Flexible parsing: scans ALL sheets for Sample Name, Target Name, Ct columns
  - Handles any Excel format with cell-by-cell scanning
  - User selects control gene and reference group after data scan
  - Calculates relative gene expression (2^-ΔΔCt)
  - Statistical analysis and downloadable results

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
# FLEXIBLE DATA SCANNING FUNCTIONS
# ------------------------------------------------------------

def scan_for_columns(df):
    """
    Scan DataFrame for Sample Name, Target Name, and Ct columns.
    Returns dict with column names if found, None otherwise.
    """
    
    # Define patterns for each required column
    patterns = {
        'sample': ['sample name', 'sample', 'samplename', 'name'],
        'target': ['target name', 'target', 'targetname', 'gene', 'assay', 'reporter'],
        'ct': ['ct', 'cт', 'cq', 'c t', 'c т']  # Handles Latin ct, Cyrillic т, and cq
    }
    
    found_columns = {}
    
    # Convert all column names to lowercase for matching
    cols_lower = {col: str(col).lower().strip().replace('_', ' ') for col in df.columns}
    
    # Find Sample Name column
    for col, col_lower in cols_lower.items():
        for pattern in patterns['sample']:
            if pattern in col_lower:
                found_columns['sample_name'] = col
                break
        if 'sample_name' in found_columns:
            break
    
    # Find Target Name column
    for col, col_lower in cols_lower.items():
        for pattern in patterns['target']:
            if pattern in col_lower:
                found_columns['target_name'] = col
                break
        if 'target_name' in found_columns:
            break
    
    # Find Ct column (special handling for exact match due to short name)
    for col, col_lower in cols_lower.items():
        # Check for exact matches or very close matches
        col_clean = col_lower.replace(' ', '').replace('_', '')
        for pattern in patterns['ct']:
            pattern_clean = pattern.replace(' ', '')
            if col_clean == pattern_clean or pattern_clean in col_clean:
                found_columns['ct'] = col
                break
        if 'ct' in found_columns:
            break
    
    # Return found columns if we have at least Sample Name, Target Name, and Ct
    if all(key in found_columns for key in ['sample_name', 'target_name', 'ct']):
        return found_columns
    else:
        return None


def parse_ct_value(value):
    """Convert Ct value to numeric, handle Undetermined"""
    if pd.isna(value):
        return np.nan
    
    s = str(value).strip().lower()
    
    # Handle undetermined values
    if s in ['undetermined', 'undet', 'na', 'n/a', '-', '', 'undeter']:
        return np.nan
    
    # Try to convert to float
    try:
        # Remove non-numeric characters except decimal point and minus
        import re
        s_clean = re.sub(r'[^0-9\.\-]', '', s)
        return float(s_clean)
    except (ValueError, TypeError):
        return np.nan


def scan_excel_file(uploaded_file):
    """
    Scan all sheets in Excel file and extract data from columns:
    Sample Name, Target Name, Ct
    
    Returns DataFrame with standardized columns.
    """
    file_name = uploaded_file.name
    all_data = []
    
    try:
        # Read all sheets
        if file_name.endswith('.csv'):
            sheets = {'Sheet1': pd.read_csv(uploaded_file)}
        else:
            sheets = pd.read_excel(uploaded_file, sheet_name=None, engine='openpyxl')
        
        for sheet_name, df in sheets.items():
            if df.empty:
                continue
            
            # Try to scan for required columns
            found_cols = scan_for_columns(df)
            
            if found_cols is None:
                # Try using first row as header if not found
                if len(df) > 0:
                    # Set first row as header and try again
                    new_header = df.iloc[0]
                    df = df[1:]
                    df.columns = new_header
                    df = df.reset_index(drop=True)
                    
                    found_cols = scan_for_columns(df)
            
            if found_cols is None:
                continue  # Skip this sheet
            
            # Extract the three required columns
            extracted_data = pd.DataFrame()
            extracted_data['Sample Name'] = df[found_cols['sample_name']].astype(str).str.strip()
            extracted_data['Target Name'] = df[found_cols['target_name']].astype(str).str.strip().str.upper()
            extracted_data['Ct'] = df[found_cols['ct']].apply(parse_ct_value)
            
            # Add metadata
            extracted_data['File'] = file_name
            extracted_data['Sheet'] = sheet_name
            
            # Remove rows with missing critical data
            extracted_data = extracted_data.dropna(subset=['Ct'])
            extracted_data = extracted_data[extracted_data['Sample Name'].str.len() > 0]
            extracted_data = extracted_data[extracted_data['Target Name'].str.len() > 0]
            extracted_data = extracted_data[extracted_data['Sample Name'] != 'nan']
            extracted_data = extracted_data[extracted_data['Target Name'] != 'NAN']
            
            if not extracted_data.empty:
                all_data.append(extracted_data)
    
    except Exception as e:
        st.error(f"Error reading {file_name}: {str(e)}")
        return pd.DataFrame()
    
    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame()


# ------------------------------------------------------------
# CALCULATION ENGINE
# ------------------------------------------------------------

class QPCRCalculator:
    """ΔΔCt Calculator for Relative Gene Expression"""
    
    def __init__(self, data, control_gene, reference_samples):
        """
        Args:
            data: DataFrame with [Sample Name, Target Name, Ct]
            control_gene: Housekeeping gene for normalization
            reference_samples: List of sample names for reference group
        """
        self.data = data.copy()
        self.control_gene = control_gene.upper()
        self.reference_samples = reference_samples
        self._validate_data()
    
    def _validate_data(self):
        """Validate inputs"""
        # Check control gene exists
        genes = self.data['Target Name'].unique()
        if self.control_gene not in genes:
            raise ValueError(f"Control gene '{self.control_gene}' not found. Available: {list(genes)}")
        
        # Check reference samples exist
        samples = self.data['Sample Name'].unique()
        missing = [s for s in self.reference_samples if s not in samples]
        if missing:
            raise ValueError(f"Reference samples not found: {missing}")
    
    def calculate_all(self):
        """Run complete analysis"""
        
        # Step 1: Calculate technical replicate statistics
        st.write("### 📊 Calculating Replicate Statistics...")
        
        grouped = self.data.groupby(['Sample Name', 'Target Name'])['Ct'].agg([
            ('Ct_Mean', 'mean'),
            ('Ct_SD', 'std'),
            ('Ct_Count', 'count')
        ]).reset_index()
        
        grouped['CV_Percent'] = (grouped['Ct_SD'] / grouped['Ct_Mean']) * 100
        grouped['High_Variation'] = grouped['CV_Percent'] > 5.0
        
        self.replicate_stats = grouped
        
        # Show high variation warning
        high_var = grouped[grouped['High_Variation']]
        if not high_var.empty:
            st.warning(f"⚠️ {len(high_var)} measurements with CV > 5%")
            with st.expander("View high variation samples"):
                st.dataframe(high_var[['Sample Name', 'Target Name', 'Ct_Mean', 'CV_Percent']], 
                           use_container_width=True)
        
        # Step 2: Calculate ΔCt (Target - Control Gene)
        st.write("### 📐 Calculating ΔCt Values...")
        
        # Get control gene Ct values
        control_ct = grouped[grouped['Target Name'] == self.control_gene][
            ['Sample Name', 'Ct_Mean', 'Ct_SD']
        ].rename(columns={'Ct_Mean': 'Control_Ct', 'Ct_SD': 'Control_SD'})
        
        # Get target genes (excluding control gene)
        targets = grouped[grouped['Target Name'] != self.control_gene].copy()
        
        # Merge and calculate ΔCt
        merged = targets.merge(control_ct, on='Sample Name', how='left')
        merged['Delta_Ct'] = merged['Ct_Mean'] - merged['Control_Ct']
        merged['Delta_Ct_SD'] = np.sqrt(merged['Ct_SD']**2 + merged['Control_SD']**2)
        
        self.delta_ct_data = merged
        
        # Step 3: Calculate ΔΔCt and Relative Expression
        st.write("### 🔬 Calculating Relative Expression...")
        
        # Calculate average ΔCt for reference samples
        reference_delta_ct = merged[
            merged['Sample Name'].isin(self.reference_samples)
        ].groupby('Target Name')['Delta_Ct'].mean().reset_index()
        reference_delta_ct.columns = ['Target Name', 'Reference_Delta_Ct']
        
        # Merge and calculate ΔΔCt
        merged = merged.merge(reference_delta_ct, on='Target Name', how='left')
        merged['Delta_Delta_Ct'] = merged['Delta_Ct'] - merged['Reference_Delta_Ct']
        merged['Relative_Expression'] = 2 ** (-merged['Delta_Delta_Ct'])
        merged['Relative_Expression_SD'] = merged['Relative_Expression'] * np.log(2) * merged['Delta_Ct_SD']
        
        self.final_results = merged
        
        # Step 4: Statistical Analysis
        st.write("### 📈 Performing Statistical Analysis...")
        
        stats_results = []
        
        for target in merged['Target Name'].unique():
            target_data = merged[merged['Target Name'] == target]
            
            # Get reference group ΔCt values
            ref_values = target_data[
                target_data['Sample Name'].isin(self.reference_samples)
            ]['Delta_Ct'].values
            
            for sample in target_data['Sample Name'].unique():
                # Check if this is a reference sample
                if sample in self.reference_samples:
                    stats_results.append({
                        'Target Name': target,
                        'Sample Name': sample,
                        'P_Value': 1.0,
                        'Significance': 'Reference'
                    })
                    continue
                
                # Get sample ΔCt values
                sample_values = target_data[
                    target_data['Sample Name'] == sample
                ]['Delta_Ct'].values
                
                # Perform t-test
                if len(sample_values) >= 2 and len(ref_values) >= 2:
                    try:
                        _, p_value = stats.ttest_ind(sample_values, ref_values)
                    except:
                        p_value = np.nan
                else:
                    p_value = np.nan
                
                # Determine significance
                if pd.isna(p_value):
                    sig = 'N/A'
                elif p_value < 0.001:
                    sig = '***'
                elif p_value < 0.01:
                    sig = '**'
                elif p_value < 0.05:
                    sig = '*'
                else:
                    sig = 'ns'
                
                stats_results.append({
                    'Target Name': target,
                    'Sample Name': sample,
                    'P_Value': p_value,
                    'Significance': sig
                })
        
        self.statistics = pd.DataFrame(stats_results)
        
        # Merge statistics
        self.final_results = self.final_results.merge(
            self.statistics,
            on=['Target Name', 'Sample Name'],
            how='left'
        )
        
        return self._get_summary()
    
    def _get_summary(self):
        """Generate summary table"""
        summary = self.final_results[[
            'Sample Name', 'Target Name', 'Ct_Mean', 'Ct_SD',
            'Delta_Ct', 'Relative_Expression', 'Relative_Expression_SD',
            'P_Value', 'Significance'
        ]].copy()
        
        # Round numeric columns
        for col in ['Ct_Mean', 'Ct_SD', 'Delta_Ct', 'Relative_Expression', 
                    'Relative_Expression_SD', 'P_Value']:
            if col in summary.columns:
                summary[col] = summary[col].round(3)
        
        return summary


# ------------------------------------------------------------
# VISUALIZATION
# ------------------------------------------------------------

def create_bar_graph(data, gene):
    """Create expression bar chart"""
    gene_data = data[data['Target Name'] == gene].copy()
    
    if gene_data.empty:
        return None
    
    gene_data = gene_data.sort_values('Sample Name')
    
    samples = gene_data['Sample Name'].values
    expression = gene_data['Relative_Expression'].values
    errors = gene_data['Relative_Expression_SD'].values
    significance = gene_data['Significance'].values
    
    # Assign colors
    colors = ['#95a5a6' if sig == 'Reference' else '#3498db' for sig in significance]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x_pos = np.arange(len(samples))
    ax.bar(x_pos, expression, color=colors, alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.errorbar(x_pos, expression, yerr=errors, fmt='none', ecolor='black',
                capsize=5, capthick=1.5, linewidth=1.5)
    
    # Add significance markers
    for i, (exp, err, sig) in enumerate(zip(expression, errors, significance)):
        if sig in ['*', '**', '***']:
            y_pos = exp + err + max(expression) * 0.05
            ax.text(i, y_pos, sig, ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Relative Gene Expression (Fold Change)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_title(f'{gene} Expression', fontsize=14, fontweight='bold', pad=20)
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_ylim(bottom=0, top=max(expression + errors) * 1.2)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    return fig


# ------------------------------------------------------------
# STREAMLIT APP
# ------------------------------------------------------------

def main():
    
    st.title("🧬 qPCR Relative Gene Expression Analysis")
    st.markdown("""
    **Flexible qPCR data analyzer** - Works with any Excel format!
    - Automatically scans all sheets for required columns
    - Extracts: Sample Name, Target Name, Ct values
    - Calculates relative gene expression (2^-ΔΔCt method)
    - Statistical analysis and publication-ready graphs
    """)
    
    # Sidebar
    st.sidebar.header("📋 How It Works")
    st.sidebar.markdown("""
    ### Process:
    1. **Upload** any Excel/CSV file
    2. App **scans** all sheets for:
       - Sample Name
       - Target Name  
       - Ct (or CT, Cт, Cq)
    3. **Select** control gene & reference samples
    4. **Analyze** & download results
    
    ### Notes:
    - Handles any Excel format
    - Scans ALL sheets automatically
    - Flexible column name detection
    """)
    
    # ========================================
    # STEP 1: Upload and Scan Files
    # ========================================
    st.header("📤 Step 1: Upload Files")
    
    uploaded_files = st.file_uploader(
        "Upload qPCR data files (Excel or CSV)",
        type=['xlsx', 'xls', 'csv'],
        accept_multiple_files=True,
        help="Upload one or more files. App will scan ALL sheets automatically."
    )
    
    if uploaded_files:
        st.success(f"✅ {len(uploaded_files)} file(s) uploaded")
        
        if st.button("🔍 Scan Files for Data", type="primary"):
            with st.spinner("Scanning all sheets for Sample Name, Target Name, and Ct columns..."):
                
                all_data = []
                scan_summary = []
                
                for file in uploaded_files:
                    scanned = scan_excel_file(file)
                    
                    if not scanned.empty:
                        all_data.append(scanned)
                        scan_summary.append({
                            'File': file.name,
                            'Rows Found': len(scanned),
                            'Sheets': scanned['Sheet'].nunique()
                        })
                    else:
                        st.warning(f"⚠️ No valid data found in {file.name}")
                
                if all_data:
                    combined_data = pd.concat(all_data, ignore_index=True)
                    st.session_state['scanned_data'] = combined_data
                    
                    st.success(f"✅ Scan complete! Found {len(combined_data)} data rows")
                    
                    # Show scan summary
                    st.subheader("📊 Scan Summary")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Rows", len(combined_data))
                    with col2:
                        st.metric("Unique Samples", combined_data['Sample Name'].nunique())
                    with col3:
                        st.metric("Unique Genes", combined_data['Target Name'].nunique())
                    
                    # Show file summary
                    with st.expander("📁 Files Scanned"):
                        st.dataframe(pd.DataFrame(scan_summary), use_container_width=True)
                    
                    # Preview data
                    st.subheader("🔍 Data Preview")
                    st.dataframe(combined_data[['Sample Name', 'Target Name', 'Ct', 'File', 'Sheet']].head(20),
                               use_container_width=True)
                    
                    # Show unique values
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("**🧬 All Genes Found:**")
                        genes = sorted(combined_data['Target Name'].unique())
                        st.write(", ".join(genes))
                    
                    with col2:
                        st.markdown("**🔬 All Samples Found:**")
                        samples = sorted(combined_data['Sample Name'].unique())
                        st.write(", ".join(samples))
                else:
                    st.error("❌ No valid data found in any uploaded files")
    
    # ========================================
    # STEP 2: Select Control Gene and Reference Samples
    # ========================================
    if 'scanned_data' in st.session_state:
        st.markdown("---")
        st.header("⚙️ Step 2: Configure Analysis")
        
        scanned_data = st.session_state['scanned_data']
        available_genes = sorted(scanned_data['Target Name'].unique())
        available_samples = sorted(scanned_data['Sample Name'].unique())
        
        col1, col2 = st.columns(2)
        
        # Control Gene Selection
        with col1:
            st.subheader("🧬 Control Gene")
            st.markdown("Select the housekeeping/reference gene for normalization")
            
            control_gene = st.selectbox(
                "Control Gene (Housekeeping)",
                options=available_genes,
                help="Common: ACTB (β-actin), GAPDH, 18S, B2M"
            )
            
            if control_gene:
                gene_count = len(scanned_data[scanned_data['Target Name'] == control_gene])
                st.info(f"✅ **{control_gene}** selected ({gene_count} measurements)")
        
        # Reference Sample Selection
        with col2:
            st.subheader("📊 Reference Samples")
            st.markdown("Select samples for the control/reference group")
            
            reference_samples = st.multiselect(
                "Reference/Control Samples",
                options=available_samples,
                help="Select all samples that belong to your control group"
            )
            
            if reference_samples:
                ref_count = len(scanned_data[scanned_data['Sample Name'].isin(reference_samples)])
                st.info(f"✅ {len(reference_samples)} sample(s) selected ({ref_count} measurements)")
                
                with st.expander("View selected reference samples"):
                    for sample in reference_samples:
                        st.write(f"• {sample}")
        
        # Save configuration
        if control_gene and reference_samples:
            st.session_state['control_gene'] = control_gene
            st.session_state['reference_samples'] = reference_samples
            
            # Show what will be compared
            st.markdown("---")
            st.info(f"""
            **Analysis Setup:**
            - Control Gene: **{control_gene}**
            - Reference Group: **{len(reference_samples)} samples**
            - Test Groups: **{len(available_samples) - len(reference_samples)} samples**
            """)
    
    # ========================================
    # STEP 3: Run Analysis
    # ========================================
    if all(key in st.session_state for key in ['scanned_data', 'control_gene', 'reference_samples']):
        st.markdown("---")
        st.header("🚀 Step 3: Run Analysis")
        
        if st.button("🧮 Calculate Relative Expression", type="primary"):
            
            scanned_data = st.session_state['scanned_data']
            control_gene = st.session_state['control_gene']
            reference_samples = st.session_state['reference_samples']
            
            try:
                with st.spinner("Running ΔΔCt analysis..."):
                    
                    calculator = QPCRCalculator(
                        data=scanned_data,
                        control_gene=control_gene,
                        reference_samples=reference_samples
                    )
                    
                    summary = calculator.calculate_all()
                    
                    st.session_state['calculator'] = calculator
                    st.session_state['summary'] = summary
                    
                    st.success("✅ Analysis Complete!")
                    
                    # Results metrics
                    st.subheader("📊 Results Summary")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        sig_count = len(summary[summary['Significance'].isin(['*', '**', '***'])])
                        st.metric("Significant Results", sig_count)
                    with col2:
                        st.metric("Genes Analyzed", summary['Target Name'].nunique())
                    with col3:
                        st.metric("Samples Analyzed", summary['Sample Name'].nunique())
                    
                    # Summary table
                    st.dataframe(summary, use_container_width=True)
                    
                    # Significance breakdown
                    st.markdown("**📈 Statistical Significance:**")
                    sig_counts = summary['Significance'].value_counts()
                    for sig, count in sig_counts.items():
                        st.write(f"  • {sig}: {count}")
            
            except Exception as e:
                st.error(f"❌ Analysis failed: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
    
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
                        
                        zip_buffer = io.BytesIO()
                        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                            
                            for gene in selected_genes:
                                fig = create_bar_graph(summary, gene)
                                if fig:
                                    img_buffer = io.BytesIO()
                                    fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                                    img_buffer.seek(0)
                                    
                                    zip_file.writestr(f"{gene}_expression.png", img_buffer.getvalue())
                                    
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
                    st.warning("Please select at least one gene")


if __name__ == "__main__":
    main()