"""
qpcr_analysis_app.py — Complete qPCR Analysis Web Application
------------------------------------------------------------
Features:
  - Flexible parsing: scans ALL sheets for Sample Name, Target Name, Ct columns
  - Handles any Excel format with cell-by-cell scanning
  - User selects Ct threshold, control gene, and reference group
  - Calculates relative gene expression (2^-ΔΔCt) with statistical rigor
  - Enhanced statistical analysis and downloadable results (Excel/ZIP)

Usage:
    streamlit run qpcr_analysis_app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import io
import re
from datetime import datetime
import zipfile
from typing import Dict, Any, List, Optional, Union

# --- CONSTANTS ---
CT_CV_THRESHOLD: float = 5.0
DEFAULT_CT_CUTOFF: float = 35.0
COLUMN_PATTERNS: Dict[str, List[str]] = {
    'sample': ['sample name', 'sample', 'samplename', 'name'],
    'target': ['target name', 'target', 'targetname', 'gene', 'assay', 'reporter'],
    'ct': ['ct', 'cт', 'cq', 'c t', 'c т', 'value']
}

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

def scan_for_columns(df: pd.DataFrame) -> Optional[Dict[str, Any]]:
    """
    Scan DataFrame for Sample Name, Target Name, and Ct columns based on patterns.
    Returns dict with column names if found, None otherwise.
    """
    
    found_columns: Dict[str, str] = {}
    cols_lower = {col: str(col).lower().strip().replace('_', ' ') for col in df.columns}
    
    # Check for required columns
    for key, patterns in COLUMN_PATTERNS.items():
        for col, col_lower in cols_lower.items():
            # Check for exact matches or partial matches
            col_clean = re.sub(r'[^a-z0-9]', '', col_lower)
            for pattern in patterns:
                pattern_clean = re.sub(r'[^a-z0-9]', '', pattern)
                
                # Prioritize exact word match, then substring match
                if pattern_clean == col_clean or pattern_clean in col_clean:
                    found_columns[key] = col
                    break
            if key in found_columns:
                break

    # Return found columns if all three critical keys are present
    if all(key in found_columns for key in ['sample', 'target', 'ct']):
        return {
            'sample_name': found_columns['sample'],
            'target_name': found_columns['target'],
            'ct': found_columns['ct']
        }
    else:
        return None

def parse_ct_value(value: Any, cutoff: float) -> Union[float, np.nan]:
    """Convert Ct value to numeric, handle Undetermined and apply cutoff."""
    if pd.isna(value):
        return np.nan
    
    s = str(value).strip().lower()
    
    # 1. Handle explicit undetermined values
    if s in ['undetermined', 'undet', 'na', 'n/a', '-', '', 'undeter', 'fail', 'no amplification']:
        return np.nan
    
    # 2. Try to convert to float
    try:
        # Remove non-numeric characters except decimal point and minus
        s_clean = re.sub(r'[^0-9\.\-]', '', s)
        ct_val = float(s_clean)
        
        # 3. Apply user-defined cutoff
        if ct_val > cutoff:
            return np.nan
        
        return ct_val
        
    except (ValueError, TypeError):
        return np.nan

def scan_excel_file(uploaded_file: io.BytesIO, ct_cutoff: float) -> pd.DataFrame:
    """
    Scan all sheets in Excel file and extract data from columns.
    Returns DataFrame with standardized columns.
    """
    file_name = uploaded_file.name
    all_data: List[pd.DataFrame] = []
    
    try:
        # Read all sheets
        if file_name.endswith('.csv'):
            sheets = {'Sheet1': pd.read_csv(uploaded_file)}
        else:
            sheets = pd.read_excel(uploaded_file, sheet_name=None, engine='openpyxl')
        
        for sheet_name, df in sheets.items():
            if df.empty:
                continue
            
            original_df = df.copy()
            found_cols = scan_for_columns(df)
            
            # Fallback: Try using first row as header if no columns found
            if found_cols is None and len(df) > 0:
                new_header = original_df.iloc[0]
                df_with_new_header = original_df[1:].copy()
                df_with_new_header.columns = new_header
                df_with_new_header = df_with_new_header.reset_index(drop=True)
                
                found_cols = scan_for_columns(df_with_new_header)
                if found_cols:
                    df = df_with_new_header # Use the new DataFrame
            
            if found_cols is None:
                continue  # Skip this sheet
            
            # Extract and standardize the three required columns
            extracted_data = pd.DataFrame()
            extracted_data['Sample Name'] = df[found_cols['sample_name']].astype(str).str.strip()
            extracted_data['Target Name'] = df[found_cols['target_name']].astype(str).str.strip().str.upper()
            
            # Apply parse_ct_value with the user-defined cutoff
            extracted_data['Ct'] = df[found_cols['ct']].apply(lambda x: parse_ct_value(x, cutoff=ct_cutoff))
            
            # Add metadata
            extracted_data['File'] = file_name
            extracted_data['Sheet'] = sheet_name
            
            # Clean up: Remove rows with missing critical data
            extracted_data = extracted_data.dropna(subset=['Ct'])
            extracted_data = extracted_data[extracted_data['Sample Name'].str.len() > 0]
            extracted_data = extracted_data[extracted_data['Target Name'].str.len() > 0]
            extracted_data = extracted_data[extracted_data['Sample Name'] != 'nan']
            extracted_data = extracted_data[extracted_data['Target Name'] != 'NAN']
            
            if not extracted_data.empty:
                all_data.append(extracted_data)
    
    except Exception as e:
        st.error(f"Error reading {file_name} on sheet {sheet_name}: {str(e)}")
        return pd.DataFrame()
    
    if all_data:
        # Reset index and ensure all columns are correct types
        combined_data = pd.concat(all_data, ignore_index=True)
        combined_data['Ct'] = combined_data['Ct'].astype(float)
        return combined_data
    else:
        return pd.DataFrame()


# ------------------------------------------------------------
# CALCULATION ENGINE (Enhanced for statistical rigor)
# ------------------------------------------------------------

class QPCRCalculator:
    """ΔΔCt Calculator for Relative Gene Expression with statistical rigor"""
    
    def __init__(self, data: pd.DataFrame, control_gene: str, reference_samples: List[str]):
        self.data = data.copy()
        self.control_gene = control_gene.upper()
        self.reference_samples = reference_samples
        self.data['Is_Reference'] = self.data['Sample Name'].isin(self.reference_samples)
        
        self.replicate_stats: pd.DataFrame
        self.delta_ct_raw: pd.DataFrame
        self.delta_delta_ct_data: pd.DataFrame
        self.final_results: pd.DataFrame
    
    def calculate_all(self) -> pd.DataFrame:
        """Run complete analysis"""
        
        # Step 1: Calculate technical replicate statistics (for QC)
        st.write("### 📊 Step 1: Calculating Technical Replicate Statistics...")
        self.replicate_stats = self._calculate_replicate_stats()
        
        # Step 2: Calculate individual raw ΔCt values
        st.write("### 📐 Step 2: Calculating Raw ΔCt Values...")
        self.delta_ct_raw = self._calculate_delta_ct_raw()
        
        # Step 3: Calculate mean ΔΔCt and Relative Expression
        st.write("### 🔬 Step 3: Calculating ΔΔCt and Relative Expression...")
        self.delta_delta_ct_data = self._calculate_delta_delta_ct()
        
        # Step 4: Statistical Analysis (T-test on raw ΔCt values)
        st.write("### 📈 Step 4: Performing Statistical Analysis (T-test on raw ΔCt)...")
        self.final_results = self._perform_statistical_analysis()
        
        return self._get_summary()
        
    def _calculate_replicate_stats(self) -> pd.DataFrame:
        """Calculate technical replicate statistics (Mean, SD, CV) per Sample/Target."""
        grouped = self.data.groupby(['Sample Name', 'Target Name'])['Ct'].agg([
            ('Ct_Mean', 'mean'),
            ('Ct_SD', 'std'),
            ('Ct_Count', 'count')
        ]).reset_index()
        
        # Calculate CV
        grouped['CV_Percent'] = (grouped['Ct_SD'] / grouped['Ct_Mean']) * 100
        grouped['High_Variation'] = grouped['CV_Percent'] > CT_CV_THRESHOLD
        
        # Show high variation warning
        high_var = grouped[grouped['High_Variation']]
        if not high_var.empty:
            st.warning(f"⚠️ {len(high_var)} measurements with CV > {CT_CV_THRESHOLD:.1f}% in technical replicates.")
            with st.expander("View high variation samples (CV > 5.0%)"):
                st.dataframe(high_var[['Sample Name', 'Target Name', 'Ct_Mean', 'CV_Percent', 'Ct_Count']], 
                           use_container_width=True)
        
        return grouped
    
    def _calculate_delta_ct_raw(self) -> pd.DataFrame:
        """Calculate ΔCt for every raw measurement (Ct_Target - Ct_Control_Mean)"""
        
        # 1. Get the control gene mean Ct for each unique sample (biological unit)
        control_ct_means = self.data[self.data['Target Name'] == self.control_gene].groupby(
            'Sample Name'
        )['Ct'].mean().reset_index().rename(columns={'Ct': 'Control_Ct_Mean'})
        
        # 2. Merge control mean back into the full dataset
        raw_data_with_control = self.data.merge(control_ct_means, on='Sample Name', how='left')
        
        # 3. Remove the control gene rows from the analysis set
        analysis_data = raw_data_with_control[
            raw_data_with_control['Target Name'] != self.control_gene
        ].copy()
        
        # 4. Calculate raw Delta Ct
        analysis_data['Delta_Ct_Raw'] = analysis_data['Ct'] - analysis_data['Control_Ct_Mean']
        
        # Remove any rows where control gene mean was missing
        analysis_data = analysis_data.dropna(subset=['Control_Ct_Mean', 'Delta_Ct_Raw'])
        
        return analysis_data[['Sample Name', 'Target Name', 'Delta_Ct_Raw', 'Is_Reference']]

    def _calculate_delta_delta_ct(self) -> pd.DataFrame:
        """Aggregates raw ΔCt values, calculates ΔΔCt and Relative Expression."""
        
        # 1. Aggregate raw Delta_Ct values to get mean Delta_Ct per sample (biological unit)
        delta_ct_agg = self.delta_ct_raw.groupby(['Sample Name', 'Target Name', 'Is_Reference'])[
            'Delta_Ct_Raw'
        ].agg([
            ('Delta_Ct_Mean', 'mean'),
            ('Delta_Ct_SD', 'std'),
            ('Delta_Ct_Count', 'count')
        ]).reset_index()
        
        # 2. Calculate average ΔCt for the reference group
        reference_delta_ct_mean = delta_ct_agg[
            delta_ct_agg['Is_Reference'] == True
        ].groupby('Target Name')['Delta_Ct_Mean'].mean().reset_index()
        reference_delta_ct_mean.columns = ['Target Name', 'Reference_Delta_Ct_Mean']
        
        # 3. Merge and calculate ΔΔCt
        merged = delta_ct_agg.merge(reference_delta_ct_mean, on='Target Name', how='left')
        merged['Delta_Delta_Ct'] = merged['Delta_Ct_Mean'] - merged['Reference_Delta_Ct_Mean']
        merged['Relative_Expression'] = 2 ** (-merged['Delta_Delta_Ct'])
        
        # 4. Calculate propagated Relative Expression SD (using Delta_Ct_SD)
        # Note: This is an approximation of the error propagation (2^-(DeltaCt)*ln(2)*SD_DeltaCt)
        merged['Relative_Expression_SD'] = merged['Relative_Expression'] * np.log(2) * merged['Delta_Ct_SD']
        
        return merged
        
    def _perform_statistical_analysis(self) -> pd.DataFrame:
        """Performs two-sample T-test on raw Delta_Ct distributions."""
        
        results_list: List[Dict[str, Any]] = []
        
        for target in self.delta_ct_raw['Target Name'].unique():
            target_raw_data = self.delta_ct_raw[
                self.delta_ct_raw['Target Name'] == target
            ]
            
            # Get raw Delta_Ct values for the reference group (full distribution)
            ref_raw_delta_ct = target_raw_data[
                target_raw_data['Is_Reference'] == True
            ]['Delta_Ct_Raw'].dropna().values
            
            for sample in target_raw_data['Sample Name'].unique():
                
                # Get raw Delta_Ct values for the test sample (full distribution)
                sample_raw_delta_ct = target_raw_data[
                    target_raw_data['Sample Name'] == sample
                ]['Delta_Ct_Raw'].dropna().values
                
                p_value: float = 1.0
                sig: str = 'Reference'
                
                if sample in self.reference_samples:
                    # Self-comparison of reference against reference is 1.0 (or not run)
                    p_value = np.nan
                    sig = 'Reference'
                
                # Perform T-test if both distributions have enough points (n>=2)
                elif len(sample_raw_delta_ct) >= 2 and len(ref_raw_delta_ct) >= 2:
                    try:
                        # Two-sided independent T-test (Assumes unequal variance: Welch's T-test)
                        _, p_value = stats.ttest_ind(sample_raw_delta_ct, ref_raw_delta_ct, equal_var=False)
                        
                        if p_value < 0.001:
                            sig = '***'
                        elif p_value < 0.01:
                            sig = '**'
                        elif p_value < 0.05:
                            sig = '*'
                        else:
                            sig = 'ns'
                            
                    except Exception:
                        p_value = np.nan
                        sig = 'N/A Error'
                else:
                    p_value = np.nan
                    sig = 'N/A Replicates'
                
                results_list.append({
                    'Target Name': target,
                    'Sample Name': sample,
                    'P_Value': p_value,
                    'Significance': sig
                })
        
        stats_df = pd.DataFrame(results_list)
        
        # Merge stats with final results
        final_merged = self.delta_delta_ct_data.merge(
            stats_df,
            on=['Target Name', 'Sample Name'],
            how='left'
        )
        
        # Bring in Ct_Mean and Ct_SD for the summary
        final_merged = final_merged.merge(
            self.replicate_stats[['Sample Name', 'Target Name', 'Ct_Mean', 'Ct_SD']],
            on=['Sample Name', 'Target Name'],
            how='left'
        )
        
        return final_merged
    
    def _get_summary(self) -> pd.DataFrame:
        """Generate summary table"""
        summary = self.final_results[[
            'Sample Name', 'Target Name', 'Ct_Mean', 'Ct_SD',
            'Delta_Ct_Mean', 'Delta_Ct_SD', 'Delta_Delta_Ct',
            'Relative_Expression', 'Relative_Expression_SD',
            'P_Value', 'Significance'
        ]].copy()
        
        # Rename columns for clarity
        summary.columns = [
            'Sample Name', 'Target Name', 'Ct_Mean', 'Ct_SD',
            'Delta Ct Mean (Target-Control)', 'Delta Ct SD', 'Delta Delta Ct',
            'Fold Change (2^-DDCt)', 'Fold Change SD',
            'P-Value', 'Significance'
        ]
        
        # Round numeric columns
        for col in ['Ct_Mean', 'Ct_SD', 'Delta Ct Mean (Target-Control)', 'Delta Ct SD',
                    'Delta Delta Ct', 'Fold Change (2^-DDCt)', 'Fold Change SD', 'P-Value']:
            if col in summary.columns:
                summary[col] = summary[col].round(4)
        
        return summary


# ------------------------------------------------------------
# VISUALIZATION
# ------------------------------------------------------------

def create_bar_graph(data: pd.DataFrame, gene: str, log_scale: bool) -> Optional[plt.Figure]:
    """Create expression bar chart with error bars and significance markers"""
    gene_data = data[data['Target Name'] == gene].copy()
    
    if gene_data.empty:
        return None
    
    # Sort samples for consistent plotting order
    gene_data = gene_data.sort_values('Sample Name')
    
    samples = gene_data['Sample Name'].values
    expression = gene_data['Fold Change (2^-DDCt)'].values
    errors = gene_data['Fold Change SD'].values
    significance = gene_data['Significance'].values
    
    # Assign colors: Reference samples are gray, others are blue
    colors = ['#95a5a6' if sig == 'Reference' else '#3498db' for sig in significance]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x_pos = np.arange(len(samples))
    ax.bar(x_pos, expression, color=colors, alpha=0.8, edgecolor='black', linewidth=1.2)
    
    # Error bars (Plotting logic for linear vs log scale must be different)
    if log_scale:
        # For log scale, plot standard deviation in the log space if possible, 
        # but here we plot it on the linear scale and let the axis handle the log conversion.
        ax.errorbar(x_pos, expression, yerr=errors, fmt='none', ecolor='black',
                    capsize=5, capthick=1.5, linewidth=1.5, zorder=3)
    else:
        # Linear scale error bars
        ax.errorbar(x_pos, expression, yerr=errors, fmt='none', ecolor='black',
                    capsize=5, capthick=1.5, linewidth=1.5, zorder=3)

    # Add significance markers
    max_err = max(expression + errors) if len(expression) > 0 else 1
    
    for i, (exp, err, sig) in enumerate(zip(expression, errors, significance)):
        if sig in ['*', '**', '***']:
            # Position marker slightly above the error bar
            y_pos = exp + err + max_err * 0.05
            ax.text(i, y_pos, sig, ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    # Formatting
    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Relative Gene Expression (Fold Change)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_title(f'{gene} Expression', fontsize=14, fontweight='bold', pad=20)
    ax.axhline(y=1, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Reference Baseline (1.0)')
    ax.grid(axis='y', alpha=0.3)
    ax.legend(loc='upper right')
    
    if log_scale:
        ax.set_yscale('log', base=2)
        ax.set_ylabel('Relative Gene Expression (Log2 Fold Change)', fontsize=12, fontweight='bold')
        # Adjust y-limit for log scale visualization
        current_min = min(expression[expression > 0]) * 0.5
        current_max = max(expression + errors) * 2
        ax.set_ylim(max(current_min, 0.01), current_max)
    else:
        ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    return fig


# ------------------------------------------------------------
# STREAMLIT APP
# ------------------------------------------------------------

def main():
    
    st.title("🧬 qPCR Relative Gene Expression Analysis")
    st.markdown("""
    **Flexible qPCR data analyzer** - Works with any Excel/CSV format!
    - Automatically scans all sheets for required columns
    - Extracts: Sample Name, Target Name, Ct values
    - Calculates relative gene expression ($\mathbf{2^{-\Delta\Delta Ct}}$) with robust statistics
    """)
    
    # Sidebar Configuration (Added Ct Cutoff here)
    st.sidebar.header("⚙️ Configuration")
    
    ct_cutoff = st.sidebar.slider(
        "Ct Cutoff / Undetermined Threshold",
        min_value=25.0,
        max_value=40.0,
        value=DEFAULT_CT_CUTOFF,
        step=0.5,
        help="Any Ct value above this threshold will be treated as Undetermined (NaN)."
    )
    
    st.sidebar.markdown("""---""")
    st.sidebar.header("📋 Analysis Steps")
    st.sidebar.markdown("""
    1. **Upload** Excel/CSV files.
    2. **Scan** files and preview data.
    3. **Select** Control Gene & Reference Samples.
    4. **Calculate** results.
    5. **Download** analysis files and graphs.
    """)
    
    
    # ========================================
    # STEP 1: Upload and Scan Files
    # ========================================
    st.header("📤 Step 1: Upload and Scan Files")
    
    uploaded_files = st.file_uploader(
        "Upload qPCR data files (Excel or CSV)",
        type=['xlsx', 'xls', 'csv'],
        accept_multiple_files=True,
        help="Upload one or more files. Ct values > cutoff (default 35.0) will be filtered."
    )
    
    if uploaded_files:
        st.success(f"✅ {len(uploaded_files)} file(s) uploaded (Ct cutoff: {ct_cutoff})")
        
        # Use session state to manage scanning status
        if 'scan_performed' not in st.session_state:
            st.session_state['scan_performed'] = False
            
        if st.button("🔍 Scan and Consolidate Data", type="primary"):
            st.session_state['scan_performed'] = True
            
            with st.spinner("Scanning all sheets for Sample Name, Target Name, and Ct columns..."):
                
                all_data: List[pd.DataFrame] = []
                scan_summary: List[Dict[str, Any]] = []
                
                for file in uploaded_files:
                    scanned = scan_excel_file(file, ct_cutoff)
                    
                    if not scanned.empty:
                        all_data.append(scanned)
                        scan_summary.append({
                            'File': file.name,
                            'Valid Rows Found': len(scanned),
                            'Unique Genes': scanned['Target Name'].nunique()
                        })
                    else:
                        st.warning(f"⚠️ No valid data found in {file.name}")
                
                if all_data:
                    combined_data = pd.concat(all_data, ignore_index=True)
                    st.session_state['scanned_data'] = combined_data
                    
                    st.success(f"✅ Scan complete! Found {len(combined_data)} valid data rows across {len(uploaded_files)} files.")
                    
                    # Show scan summary
                    st.subheader("📊 Scan Summary")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Valid Rows", len(combined_data))
                    with col2:
                        st.metric("Unique Samples", combined_data['Sample Name'].nunique())
                    with col3:
                        st.metric("Unique Genes", combined_data['Target Name'].nunique())
                    
                    # Show file summary
                    with st.expander("📁 Files Consolidated"):
                        st.dataframe(pd.DataFrame(scan_summary), use_container_width=True)
                    
                    # Preview data
                    st.subheader("🔍 Data Preview (First 20 Rows)")
                    st.dataframe(combined_data[['Sample Name', 'Target Name', 'Ct', 'File']].head(20),
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
                    st.error("❌ No valid data found in any uploaded files after applying Ct cutoff.")
                    del st.session_state['scanned_data'] # Clear state if scan failed
    
    # ========================================
    # STEP 2: Select Control Gene and Reference Samples
    # ========================================
    if 'scanned_data' in st.session_state:
        st.markdown("---")
        st.header("⚙️ Step 2: Configure Analysis Parameters")
        
        scanned_data = st.session_state['scanned_data']
        available_genes = sorted(scanned_data['Target Name'].unique())
        available_samples = sorted(scanned_data['Sample Name'].unique())
        
        col1, col2 = st.columns(2)
        
        # Control Gene Selection
        with col1:
            st.subheader("🧬 Control Gene")
            st.markdown("Select the housekeeping/reference gene for $C_t$ normalization.")
            
            control_gene = st.selectbox(
                "Control Gene (Housekeeping)",
                options=available_genes,
                help="Common: ACTB, GAPDH, 18S, B2M. Note: Gene names are converted to UPPERCASE."
            )
            
        # Reference Sample Selection
        with col2:
            st.subheader("📊 Reference Samples")
            st.markdown("Select all biological samples defining the control/reference group.")
            
            reference_samples = st.multiselect(
                "Reference/Control Samples",
                options=available_samples,
                help="Select all samples that will be used to calculate the $\Delta\Delta Ct$ baseline."
            )
            
        
        # Save configuration
        if control_gene and reference_samples:
            st.session_state['control_gene'] = control_gene
            st.session_state['reference_samples'] = reference_samples
            
            # Show what will be compared
            st.markdown("---")
            st.info(f"""
            **Analysis Setup Confirmed:**
            - Control Gene: **{control_gene}**
            - Reference Group: **{len(reference_samples)} samples**
            - Comparison (Test) Samples: **{len(available_samples) - len(reference_samples)} samples**
            """)
    
    # ========================================
    # STEP 3: Run Analysis
    # ========================================
    if all(key in st.session_state for key in ['scanned_data', 'control_gene', 'reference_samples']):
        st.markdown("---")
        st.header("🚀 Step 3: Run $\Delta\Delta Ct$ Analysis")
        
        if st.button("🧮 Calculate Relative Expression & Statistics", type="primary"):
            
            scanned_data = st.session_state['scanned_data']
            control_gene = st.session_state['control_gene']
            reference_samples = st.session_state['reference_samples']
            
            # Reset results state
            st.session_state['calculator'] = None
            st.session_state['summary'] = None
            
            try:
                with st.spinner("Running ΔΔCt analysis and T-tests..."):
                    
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
                    st.subheader("📊 Results Overview")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        sig_count = len(summary[summary['Significance'].isin(['*', '**', '***'])])
                        st.metric("Significant Results (*, **, ***)", sig_count)
                    with col2:
                        st.metric("Analysed Groups", summary['Target Name'].nunique() * summary['Sample Name'].nunique())
                    with col3:
                        st.metric("Calculated T-tests", len(summary) - len(summary[summary['Significance'] == 'Reference']))
                    
                    # Summary table
                    st.subheader("Final Relative Expression Summary (Fold Change)")
                    st.markdown("This table shows the mean Fold Change and P-Value calculated via T-test on individual $\Delta Ct$ values.")
                    st.dataframe(summary, use_container_width=True)
                    
                    # Significance breakdown
                    st.markdown("**📈 Statistical Significance Legend:**")
                    st.write("  • `***`: p < 0.001")
                    st.write("  • `**`: p < 0.01")
                    st.write("  • `*`: p < 0.05")
                    st.write("  • `ns`: p ≥ 0.05 (Not Significant)")
                    st.write("  • `Reference`: The sample is part of the control group.")
            
            except Exception as e:
                st.error(f"❌ Analysis failed. Check if your control gene and reference samples contain valid Ct data: {str(e)}")
                # Optional: display traceback for advanced debugging
                # import traceback
                # st.code(traceback.format_exc())
    
    # ========================================
    # STEP 4: Download Results and Graphs
    # ========================================
    if 'calculator' in st.session_state and 'summary' in st.session_state and st.session_state['calculator']:
        st.markdown("---")
        st.header("📥 Step 4: Download Results and Graphs")
        
        calculator: QPCRCalculator = st.session_state['calculator']
        summary = st.session_state['summary']
        
        col1, col2 = st.columns(2)
        
        # Excel Download
        with col1:
            st.subheader("📊 Download Analysis Files")
            
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                calculator.replicate_stats.to_excel(writer, sheet_name='1_Replicate_QC', index=False)
                calculator.delta_ct_raw.to_excel(writer, sheet_name='2_Delta_Ct_Raw', index=False)
                calculator.delta_delta_ct_data.to_excel(writer, sheet_name='3_Delta_Delta_Ct_Data', index=False)
                summary.to_excel(writer, sheet_name='4_Final_Summary', index=False)
            
            excel_data = output.getvalue()
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            st.download_button(
                label="📥 Download All Excel Results (.xlsx)",
                data=excel_data,
                file_name=f"qpcr_results_complete_{timestamp}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        # Graph Generation
        with col2:
            st.subheader("📈 Generate and Download Graphs")
            
            genes = summary['Target Name'].unique()
            selected_genes = st.multiselect(
                "Select genes to plot",
                options=genes,
                default=list(genes),
                help="Select which target genes you want to visualize."
            )
            
            log_scale_toggle = st.checkbox("Plot Y-axis on Log2 Scale (Recommended)", value=True)
            
            if st.button("Generate Graphs"):
                if selected_genes:
                    with st.spinner("Creating graphs..."):
                        
                        zip_buffer = io.BytesIO()
                        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                            
                            for gene in selected_genes:
                                fig = create_bar_graph(summary, gene, log_scale_toggle)
                                if fig:
                                    img_buffer = io.BytesIO()
                                    fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                                    img_buffer.seek(0)
                                    
                                    # Write PNG to zip file
                                    zip_file.writestr(f"{gene}_expression.png", img_buffer.getvalue())
                                    
                                    # Display in Streamlit
                                    st.pyplot(fig)
                                    plt.close(fig)
                        
                        zip_buffer.seek(0)
                        
                        # Download button for the ZIP file
                        st.download_button(
                            label="📥 Download All Graphs (ZIP)",
                            data=zip_buffer.getvalue(),
                            file_name=f"qpcr_graphs_{timestamp}.zip",
                            mime="application/zip"
                        )
                else:
                    st.warning("Please select at least one gene to plot.")


if __name__ == "__main__":
    # Ensure session state variables are initialized
    if 'scanned_data' not in st.session_state:
        st.session_state['scanned_data'] = None
    if 'calculator' not in st.session_state:
        st.session_state['calculator'] = None
    if 'summary' not in st.session_state:
        st.session_state['summary'] = None
        
    main()
