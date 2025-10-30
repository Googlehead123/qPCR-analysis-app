"""
qpcr_analysis_app.py — Complete qPCR Analysis Web Application
------------------------------------------------------------
Features:
  - Robust Dynamic Parsing: Scans all rows/sheets to find data table start.
  - NEW: Strict matching for 'Sample Name' and 'Target Name'.
  - NEW: Lenient matching for 'CT' variations (Ct, Cт, Cq, etc.).
  - Multi-File Consolidation: Handles multiple uploaded files as a single experiment.
  - Ct Cutoff: User-configurable threshold for treating data as Undetermined (NaN).
  - Relative Expression (2^-ΔΔCt) Calculation.
  - Statistically Rigorous T-test on individual Delta Ct values.
  - Downloadable results (Excel/ZIP) and graphs.

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
from typing import Dict, Any, List, Optional, Union, Tuple

# --- CONSTANTS ---
CT_CV_THRESHOLD: float = 5.0
DEFAULT_CT_CUTOFF: float = 40.0

# Define patterns for required columns and their strictness level
# 'strict': True requires the cleaned cell value (lower/stripped) to exactly match the pattern.
# 'strict': False uses broader matching (contains/similar after removing special chars).
COLUMN_CONFIG: Dict[str, Dict[str, Any]] = {
    'sample': {'patterns': ['sample name'], 'strict': True},
    'target': {'patterns': ['target name'], 'strict': True},
    # Lenient patterns for Ct, including common variations and the Cyrillic 'т'
    'ct': {'patterns': ['ct', 'cт', 'cq', 'c t', 'c т', 'value', 'ct value', 'qpcr ct'], 'strict': False},
}


# Set page config
st.set_page_config(
    page_title="qPCR Analysis Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ------------------------------------------------------------
# ROBUST DATA SCANNING FUNCTIONS (CRITICAL UPDATE)
# ------------------------------------------------------------

def scan_for_columns_and_start(df: pd.DataFrame, file_name: str) -> Tuple[Optional[Dict[str, str]], int]:
    """
    Scan all rows of a DataFrame to dynamically find the header row (start_row)
    that contains the required column names based on strictness rules.
    
    Returns: (found_columns_map, start_row_index)
    """
    
    # Try the first 100 rows to find a viable header
    for i in range(min(100, len(df))):
        
        # 1. Get the raw values of the potential header row (pandas uses numeric indices when header=None)
        header_row_raw = df.iloc[i].astype(str).tolist()
        
        found_columns: Dict[str, str] = {}
        
        # 2. Iterate through each cell in the potential header row
        for col_index, cell_value in enumerate(header_row_raw):
            # Standard cleaning: lower and strip whitespace
            cell_clean = cell_value.strip().lower()

            # 3. Check against required patterns
            for required_key, config in COLUMN_CONFIG.items():
                
                # If we already found this required column, skip to next required_key
                if required_key in found_columns:
                    continue
                
                # Check patterns
                for pattern in config['patterns']:
                    pattern_clean = pattern.strip().lower()
                    
                    if config['strict']:
                        # STRICT MATCH (e.g., 'sample name', 'target name'): Must be an exact match
                        if cell_clean == pattern_clean:
                            # Store the original column index's title (which is the numeric index if header=None)
                            # This is the actual column identifier in the raw DataFrame
                            found_columns[required_key] = df.columns[col_index]
                            break 
                    else:
                        # LENIENT MATCH (e.g., 'ct'): Checks for similarity/contains after stripping non-alphanumeric chars
                        cell_alnum = re.sub(r'[^a-z0-9]', '', cell_clean)
                        pattern_alnum = re.sub(r'[^a-z0-9]', '', pattern_clean)
                        
                        if pattern_alnum == cell_alnum or pattern_alnum in cell_alnum:
                            found_columns[required_key] = df.columns[col_index]
                            break
        
        # Check if all critical columns were found in this row
        if all(key in found_columns for key in ['sample', 'target', 'ct']):
            # SUCCESS: Return the map of standard names to the raw index/title, and the row number
            return {
                'sample_name': found_columns['sample'],
                'target_name': found_columns['target'],
                'ct': found_columns['ct']
            }, i
    
    # Failure: Fallback if no clear header is found
    return None, -1

def parse_ct_value(value: Any, cutoff: float) -> Union[float, np.nan]:
    """Convert Ct value to numeric, handle Undetermined and apply cutoff."""
    if pd.isna(value):
        return np.nan
    
    s = str(value).strip().lower()
    
    # 1. Handle explicit undetermined values
    if s in ['undetermined', 'undet', 'na', 'n/a', '-', '', 'undeter', 'fail', 'no amplification', 'n o a m p', 'n.a.']:
        return np.nan
    
    # 2. Try to convert to float
    try:
        # Remove non-numeric characters except decimal point and minus
        s_clean = re.sub(r'[^0-9\.\-]', '', s)
        ct_val = float(s_clean)
        
        # 3. Apply user-defined cutoff
        if ct_val >= cutoff: # Use >= for threshold
            return np.nan
        
        return ct_val
        
    except (ValueError, TypeError):
        return np.nan

def scan_excel_file(uploaded_file: io.BytesIO, ct_cutoff: float) -> pd.DataFrame:
    """
    Reads an Excel/CSV file, finds the header row, and extracts required data.
    """
    file_name = uploaded_file.name
    all_sheet_data: List[pd.DataFrame] = []
    
    try:
        # Read all sheets/file content without an initial header to scan all rows
        if file_name.lower().endswith('.csv'):
             # Use the more robust C engine for CSV unless a specific delimiter is needed, 
             # but keep header=None for full scanning
            sheets = {'Sheet1': pd.read_csv(uploaded_file, header=None, encoding='utf-8', on_bad_lines='skip')}
        else:
            sheets = pd.read_excel(uploaded_file, sheet_name=None, header=None, engine='openpyxl')
        
        for sheet_name, df_raw in sheets.items():
            if df_raw.empty:
                continue
            
            # --- Dynamic Header Detection ---
            found_cols_map, header_row_index = scan_for_columns_and_start(df_raw, file_name)
            
            if found_cols_map is None or header_row_index == -1:
                st.info(f"Skipping sheet '{sheet_name}' in **{file_name}**: Could not find the required header row (Sample Name, Target Name, and Ct variation) within the first 100 rows.")
                continue
            
            # --- Extract Data using the detected column indices ---
            
            # Get the column index names (which are just numbers 0, 1, 2... because we used header=None)
            sample_col_index = found_cols_map['sample_name']
            target_col_index = found_cols_map['target_name']
            ct_col_index = found_cols_map['ct']

            # Slice the data: Rows below the header row
            df_data = df_raw[header_row_index + 1:].copy()
            
            extracted_data = pd.DataFrame()
            
            # 1. Sample Name (String, clean whitespace) - using the numeric index found
            extracted_data['Sample Name'] = df_data[sample_col_index].astype(str).str.strip()
            
            # 2. Target Name (String, uppercase and clean whitespace) - using the numeric index found
            extracted_data['Target Name'] = df_data[target_col_index].astype(str).str.strip().str.upper()
            
            # 3. Ct (Numeric, apply parsing and cutoff) - using the numeric index found
            extracted_data['Ct'] = df_data[ct_col_index].apply(lambda x: parse_ct_value(x, cutoff=ct_cutoff))
            
            # --- Final Cleaning and Metadata ---
            extracted_data['File'] = file_name
            extracted_data['Sheet'] = sheet_name
            
            # Remove rows with missing critical data (Sample, Target, or filtered Ct)
            extracted_data = extracted_data.dropna(subset=['Ct'])
            extracted_data = extracted_data[extracted_data['Sample Name'].str.len() > 0]
            extracted_data = extracted_data[extracted_data['Target Name'].str.len() > 0]
            
            if not extracted_data.empty:
                all_sheet_data.append(extracted_data)
    
    except Exception as e:
        # Provide a more generic error if parsing fails unexpectedly
        st.error(f"Error reading/parsing **{file_name}** on sheet **{sheet_name}**. The file might be corrupted or in an incompatible format. Details: {str(e)}")
        return pd.DataFrame()
    
    if all_sheet_data:
        # Concatenate data from all valid sheets/files
        combined_data = pd.concat(all_sheet_data, ignore_index=True)
        # Ensure final data types are correct
        combined_data['Ct'] = combined_data['Ct'].astype(float)
        return combined_data
    else:
        return pd.DataFrame()


# ------------------------------------------------------------
# CALCULATION ENGINE (No changes needed here)
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
# VISUALIZATION (No changes needed here)
# ------------------------------------------------------------

def create_bar_graph(data: pd.DataFrame, gene: str, log_scale: bool) -> Optional[plt.Figure]:
    """Create expression bar chart with error bars and significance markers"""
    gene_data = data[data['Target Name'] == gene].copy()
    
    if gene_data.empty:
        return None
    
    # Sort samples for consistent plotting order
    ref_samples = gene_data[gene_data['Significance'] == 'Reference']['Sample Name'].unique()
    other_samples = gene_data[gene_data['Significance'] != 'Reference']['Sample Name'].unique()
    
    sample_order = list(ref_samples) + list(other_samples)
    
    # Reindex the DataFrame according to the custom order
    gene_data['Sample_Order'] = pd.Categorical(gene_data['Sample Name'], categories=sample_order, ordered=True)
    gene_data = gene_data.sort_values('Sample_Order')
    
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
    
    # Error bars 
    ax.errorbar(x_pos, expression, yerr=errors, fmt='none', ecolor='black',
                capsize=5, capthick=1.5, linewidth=1.5, zorder=3)

    # Add significance markers
    max_err_heights = (expression + errors).replace([np.inf, -np.inf], np.nan).dropna()
    max_plot_height = max(max_err_heights) if len(max_err_heights) > 0 else 1
    
    for i, (exp, err, sig) in enumerate(zip(expression, errors, significance)):
        if sig in ['*', '**', '***']:
            y_pos = exp + err + max_plot_height * 0.05
            ax.text(i, y_pos, sig, ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    # Formatting
    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Relative Gene Expression (Fold Change)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_title(f'{gene} Expression (Normalized to {st.session_state.get("control_gene", "Control")})', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.axhline(y=1, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Reference Baseline (1.0)')
    ax.grid(axis='y', alpha=0.3)
    ax.legend(loc='upper right')
    
    if log_scale:
        ax.set_yscale('log', base=2)
        ax.set_ylabel('Relative Gene Expression (Log₂ Fold Change)', fontsize=12, fontweight='bold')
        valid_expression = expression[expression > 0]
        current_min = min(valid_expression) * 0.5 if len(valid_expression) > 0 else 0.5
        current_max = max(expression + errors) * 2
        ax.set_ylim(max(current_min, 0.01), current_max)
    else:
        ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    return fig


# ------------------------------------------------------------
# STREAMLIT APP (No changes needed here)
# ------------------------------------------------------------

def main():
    
    st.title("🧬 Robust qPCR Relative Gene Expression Analysis")
    st.markdown("""
    **Flexible qPCR data analyzer** - Upload single or multiple files. The app automatically finds the data table,
    consolidates results, and performs $\mathbf{2^{-\Delta\Delta Ct}}$ analysis with rigorous statistics.
    """)
    
    # Sidebar Configuration (Added Ct Cutoff here)
    st.sidebar.header("⚙️ Configuration")
    
    ct_cutoff = st.sidebar.slider(
        "Ct Cutoff / Undetermined Threshold",
        min_value=25.0,
        max_value=40.0,
        value=DEFAULT_CT_CUTOFF,
        step=0.5,
        help="Any Ct value **equal to or above** this threshold will be treated as Undetermined (NaN) and filtered out."
    )
    
    st.sidebar.markdown("""---""")
    st.sidebar.header("📋 Analysis Steps")
    st.sidebar.markdown("""
    1. **Upload** Excel/CSV files.
    2. **Scan** files and verify consolidated data.
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
        help="Upload one or more files for consolidated analysis."
    )
    
    if uploaded_files:
        
        # Use a unique key to force a rescan if files change
        scan_button_key = f"scan_button_{len(uploaded_files)}_{ct_cutoff}"
        
        if st.button("🔍 Scan and Consolidate Data", type="primary", key=scan_button_key):
            
            with st.spinner(f"Scanning {len(uploaded_files)} file(s) for data start rows..."):
                
                all_data: List[pd.DataFrame] = []
                scan_summary: List[Dict[str, Any]] = []
                
                for file in uploaded_files:
                    # Reset the file stream pointer for reading
                    file.seek(0)
                    scanned = scan_excel_file(file, ct_cutoff)
                    
                    if not scanned.empty:
                        all_data.append(scanned)
                        scan_summary.append({
                            'File': file.name,
                            'Valid Rows Found': len(scanned),
                            'Unique Genes': scanned['Target Name'].nunique()
                        })
                    else:
                        st.warning(f"⚠️ No valid data found in **{file.name}** after applying Ct cutoff or required columns were missing.")
                
                if all_data:
                    combined_data = pd.concat(all_data, ignore_index=True)
                    st.session_state['scanned_data'] = combined_data
                    st.session_state['analysis_configured'] = False # Reset subsequent steps
                    st.session_state['analysis_run'] = False 
                    
                    st.success(f"✅ Scan complete! Found {len(combined_data)} valid data rows across {len(all_data)} file(s).")
                    
                    # Show scan summary
                    st.subheader("📊 Consolidated Data Summary")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Valid Rows", len(combined_data))
                    with col2:
                        st.metric("Unique Samples", combined_data['Sample Name'].nunique())
                    with col3:
                        st.metric("Unique Genes", combined_data['Target Name'].nunique())
                    
                    # Show file summary
                    with st.expander("📁 Files Scanned Successfully"):
                        st.dataframe(pd.DataFrame(scan_summary), use_container_width=True)
                    
                    # Preview data for verification (user request)
                    st.subheader("🔍 Consolidated Data Preview (First 20 Rows)")
                    st.markdown("Verify the correct **Sample Name**, **Target Name**, and **Ct** values are shown.")
                    st.dataframe(combined_data[['Sample Name', 'Target Name', 'Ct', 'File']].head(20),
                               use_container_width=True)
                else:
                    st.error("❌ No valid data found in any uploaded files after applying Ct cutoff.")
                    if 'scanned_data' in st.session_state:
                        del st.session_state['scanned_data']
    
    # ========================================
    # STEP 2: Select Control Gene and Reference Samples
    # ========================================
    if 'scanned_data' in st.session_state and st.session_state['scanned_data'] is not None:
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
            
            default_gene = st.session_state.get('control_gene') if st.session_state.get('control_gene') in available_genes else (available_genes[0] if available_genes else None)
            
            control_gene = st.selectbox(
                "Control Gene (Housekeeping)",
                options=available_genes,
                index=available_genes.index(default_gene) if default_gene else 0,
                help="Common: ACTB, GAPDH, 18S, B2M. Note: Gene names are converted to UPPERCASE."
            )
            
        # Reference Sample Selection
        with col2:
            st.subheader("📊 Reference Samples")
            st.markdown("Select all biological samples defining the control/reference group.")
            
            default_samples = st.session_state.get('reference_samples', [])
            
            reference_samples = st.multiselect(
                "Reference/Control Samples",
                options=available_samples,
                default=[s for s in default_samples if s in available_samples],
                help="Select all samples that will be used to calculate the $\Delta\Delta Ct$ baseline."
            )
            
        
        # Save configuration
        if control_gene and reference_samples:
            st.session_state['control_gene'] = control_gene
            st.session_state['reference_samples'] = reference_samples
            st.session_state['analysis_configured'] = True
            st.session_state['analysis_run'] = False # Reset analysis run if config changes
            
            # Show what will be compared
            st.markdown("---")
            st.info(f"""
            **Analysis Setup Confirmed:**
            - Control Gene: **{control_gene}**
            - Reference Group: **{len(reference_samples)} samples**
            - Comparison (Test) Samples: **{len(available_samples) - len(reference_samples)} samples**
            """)
        else:
            st.warning("Please select a **Control Gene** and at least one **Reference Sample** to continue.")
            st.session_state['analysis_configured'] = False


    # ========================================
    # STEP 3: Run Analysis
    # ========================================
    if st.session_state.get('analysis_configured', False):
        st.markdown("---")
        st.header("🚀 Step 3: Run $\Delta\Delta Ct$ Analysis")
        
        if st.button("🧮 Calculate Relative Expression & Statistics", type="primary"):
            
            scanned_data = st.session_state['scanned_data']
            control_gene = st.session_state['control_gene']
            reference_samples = st.session_state['reference_samples']
            
            # Reset results state
            st.session_state['calculator'] = None
            st.session_state['summary'] = None
            st.session_state['analysis_run'] = False
            
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
                    st.session_state['analysis_run'] = True
                    
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


    # ========================================
    # STEP 4: Download Results and Graphs
    # ========================================
    if st.session_state.get('analysis_run', False):
        st.markdown("---")
        st.header("📥 Step 4: Download Results and Graphs")
        
        calculator: QPCRCalculator = st.session_state['calculator']
        summary = st.session_state['summary']
        
        col1, col2 = st.columns(2)
        
        # Excel Download
        with col1:
            st.subheader("📊 Download Analysis Files")
            
            output = io.BytesIO()
            calc_data = calculator.replicate_stats.copy()
            
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
            
            log_scale_toggle = st.checkbox("Plot Y-axis on Log₂ Scale (Recommended)", value=True)
            
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
                                else:
                                    st.warning(f"Could not generate graph for gene: {gene} (Missing data).")

                        
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
    if 'analysis_configured' not in st.session_state:
        st.session_state['analysis_configured'] = False
    if 'analysis_run' not in st.session_state:
        st.session_state['analysis_run'] = False
        
    main()
