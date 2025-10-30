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
import re

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
        'ct': ['ct', 'cт', 'cq', 'quantification cycle'],
    }
    
    col_map = {}
    
    # Normalize column names for case-insensitive matching and removal of non-alphanumeric
    normalized_cols = {re.sub(r'[^a-z0-9]', '', str(col).lower()): col for col in df.columns}
    
    # Attempt to find the best match for each required column
    for key, search_terms in patterns.items():
        found = False
        for term in search_terms:
            normalized_term = re.sub(r'[^a-z0-9]', '', term)
            for norm_col, original_col in normalized_cols.items():
                if norm_col == normalized_term or normalized_term in norm_col:
                    col_map[key] = original_col
                    found = True
                    break
            if found:
                break
    
    required_cols = ['sample', 'target', 'ct']
    if all(k in col_map for k in required_cols):
        return {k: col_map[k] for k in required_cols}
    else:
        # Report which columns are missing
        missing = [k for k in required_cols if k not in col_map]
        st.error(f"Could not find required columns. Missing: {', '.join(missing)}. Check the data format.")
        return None

def parse_qpcr_data(uploaded_file, file_extension):
    """
    Loads data from uploaded file (Excel or CSV), finds the header, 
    and concatenates data from all sheets (if Excel).
    """
    
    try:
        if file_extension in ['xlsx', 'xls']:
            # Read all sheets from Excel file
            xls = pd.ExcelFile(uploaded_file)
            sheet_names = xls.sheet_names
            data_frames = []
            
            for sheet_name in sheet_names:
                df = pd.read_excel(xls, sheet_name=sheet_name, header=None)
                # Find header row
                header_row_index = -1
                col_names = None
                for i in range(len(df)):
                    temp_df = df.iloc[i:].copy()
                    temp_df.columns = temp_df.iloc[0]
                    temp_df = temp_df[1:].reset_index(drop=True)
                    
                    col_names_map = scan_for_columns(temp_df)
                    if col_names_map:
                        header_row_index = i
                        col_names = col_names_map
                        break
                
                if header_row_index != -1:
                    df = pd.read_excel(xls, sheet_name=sheet_name, header=header_row_index)
                    data_frames.append(df)
            
            if not data_frames:
                st.error("Could not find a valid data header (with Sample Name, Target Name, and Ct columns) in any sheet.")
                return None, None
                
            all_data = pd.concat(data_frames, ignore_index=True)
            col_names = scan_for_columns(all_data) # Re-scan columns on combined data
            
        elif file_extension == 'csv':
            df = pd.read_csv(uploaded_file, header=None)
            header_row_index = -1
            col_names = None
            
            # Find header row for CSV
            for i in range(len(df)):
                temp_df = df.iloc[i:].copy()
                temp_df.columns = temp_df.iloc[0]
                temp_df = temp_df[1:].reset_index(drop=True)
                
                col_names_map = scan_for_columns(temp_df)
                if col_names_map:
                    header_row_index = i
                    col_names = col_names_map
                    break

            if header_row_index != -1:
                all_data = pd.read_csv(uploaded_file, header=header_row_index)
                col_names = scan_for_columns(all_data)
            else:
                st.error("Could not find a valid data header (with Sample Name, Target Name, and Ct columns) in the CSV file.")
                return None, None
        else:
            st.error("Unsupported file format.")
            return None, None

        if col_names:
            return all_data, col_names
        else:
            st.error("Could not identify all necessary columns (Sample Name, Target Name, Ct) after data loading.")
            return None, None

    except Exception as e:
        st.error(f"An error occurred during file parsing: {e}")
        return None, None

# ------------------------------------------------------------
# DATA CLEANING AND STANDARDIZATION
# ------------------------------------------------------------

def standardize_sample_names(df, sample_col):
    """
    REVISION: Attempts to standardize numeric sample names (e.g., convert 1.0 to 1) 
    to ensure they are treated as the same group for analysis.
    """
    # Create a new series based on the sample column, forcing all to be strings for consistency
    s = df[sample_col].astype(str).str.strip().copy()
    
    # Attempt to convert to numeric, coercing non-convertible values to NaN
    numeric_s = pd.to_numeric(s, errors='coerce')
    
    # Check which values are non-NaN (i.e., successfully converted to a number)
    is_numeric = numeric_s.notna()
    
    # For numeric values, convert to integer (rounding first) and then back to string.
    # This turns '1.0', 1.0, 1, '1' into the string '1'.
    # We use a try/except block just in case the conversion to int fails on the rounded float (though unlikely)
    try:
        df.loc[is_numeric, sample_col] = numeric_s[is_numeric].round(0).astype(int).astype(str)
    except:
        # Fallback to just converting back to string if int conversion fails for some reason
        df.loc[is_numeric, sample_col] = numeric_s[is_numeric].astype(str).str.replace(r'\.0$', '', regex=True)
    
    # Ensure all names are stripped of whitespace
    df[sample_col] = df[sample_col].astype(str).str.strip()
    
    return df

def clean_data(df, col_names):
    """
    Filters out non-numeric Ct values, ensures columns are correct types, 
    and handles potential NaN values.
    """
    # Rename columns for simpler internal use
    df = df.rename(columns={
        col_names['sample']: 'Sample Name',
        col_names['target']: 'Target Name',
        col_names['ct']: 'Ct'
    })
    
    # Convert Ct to numeric, setting errors='coerce' to turn invalid entries into NaN
    df['Ct'] = pd.to_numeric(df['Ct'], errors='coerce')
    
    # Drop rows where Ct is NaN (i.e., failed, not detected, or invalid)
    df = df.dropna(subset=['Ct'])
    
    # Filter out Ct values that are too high (e.g., > 35) as they are typically non-specific/undetected
    df = df[df['Ct'] <= 35] 

    # Ensure Target Name and Sample Name are strings
    df['Target Name'] = df['Target Name'].astype(str).str.strip()
    df['Sample Name'] = df['Sample Name'].astype(str).str.strip()
    
    # Drop any rows where Sample Name or Target Name might be empty strings after stripping
    df = df[df['Sample Name'].str.len() > 0]
    df = df[df['Target Name'].str.len() > 0]
    
    return df[['Sample Name', 'Target Name', 'Ct']]

# ------------------------------------------------------------
# CALCULATION FUNCTIONS
# ------------------------------------------------------------

def calculate_ddct(df, hk_gene, ref_group):
    """
    Performs the 2^-ΔΔCt calculation.
    """
    
    # 1. Calculate ΔCt (Ct_Target - Ct_HK) for each replicate
    # Group by Sample and Target, calculate mean Ct for each
    ct_mean = df.groupby(['Sample Name', 'Target Name'])['Ct'].mean().reset_index()
    ct_mean = ct_mean.rename(columns={'Ct': 'Ct_Mean'})
    
    # Merge Ct_Mean back to the original filtered dataframe to get individual ΔCt calculations
    df_merged = pd.merge(df, ct_mean, on=['Sample Name', 'Target Name'])

    # Get Housekeeping (HK) mean Ct for each Sample
    hk_ct = df_merged[df_merged['Target Name'] == hk_gene].groupby('Sample Name')['Ct_Mean'].first().reset_index()
    hk_ct = hk_ct.rename(columns={'Ct_Mean': 'Ct_HK_Mean'})
    
    # Merge HK Ct back to the main dataframe
    df_ddct = pd.merge(df_merged, hk_ct, on='Sample Name', how='left')

    # Calculate ΔCt for each biological/technical replicate
    df_ddct['dCt'] = df_ddct['Ct'] - df_ddct['Ct_HK_Mean']
    
    # Calculate Mean ΔCt and SE ΔCt for each (Sample, Target) pair (n=number of replicates)
    dCt_summary = df_ddct.groupby(['Sample Name', 'Target Name'])['dCt'].agg(
        dCt_Mean='mean',
        dCt_SE=lambda x: stats.sem(x) if len(x) > 1 else 0, # Standard Error of the Mean
        N='count'
    ).reset_index()
    
    # 2. Calculate ΔΔCt (ΔCt_Sample - ΔCt_Reference)
    
    # Get the ΔCt_Ref_Mean for each target gene
    ref_dCt = dCt_summary[dCt_summary['Sample Name'] == ref_group][['Target Name', 'dCt_Mean']]
    ref_dCt = ref_dCt.rename(columns={'dCt_Mean': 'dCt_Ref_Mean'})
    
    # Merge reference ΔCt back to the summary dataframe
    summary = pd.merge(dCt_summary, ref_dCt, on='Target Name', how='left')
    
    # Calculate ΔΔCt
    summary['ddCt'] = summary['dCt_Mean'] - summary['dCt_Ref_Mean']
    
    # 3. Calculate Relative Expression (2^-ΔΔCt)
    summary['Relative Expression (2^-ddCt)'] = 2**(-summary['ddCt'])
    
    # Calculate Standard Deviation of 2^-ΔΔCt (Max and Min fold change)
    # 2^-(ΔΔCt + SE) and 2^-(ΔΔCt - SE) gives the min and max fold change
    # Note: We use dCt_SE here which is the standard error of the mean dCt. 
    # This is an approximation for error propagation in ddCt methods.
    summary['Fold Change Min'] = 2**(-(summary['ddCt'] + summary['dCt_SE']))
    summary['Fold Change Max'] = 2**(-(summary['ddCt'] - summary['dCt_SE']))
    
    # Final result includes all calculated values
    final_summary = summary.drop(columns=['dCt_Ref_Mean']).rename(columns={'dCt_Mean': 'dCt_Sample_Mean', 'dCt_SE': 'dCt_SE'})
    
    return final_summary

# ------------------------------------------------------------
# GRAPHING FUNCTIONS
# ------------------------------------------------------------

def create_bar_graph(summary_df, gene):
    """
    Creates a bar graph for a single gene across all samples.
    """
    gene_data = summary_df[summary_df['Target Name'] == gene].copy()
    
    if gene_data.empty:
        return None

    # Calculate the error bar extent (upper and lower bounds)
    gene_data['Error_Lower'] = gene_data['Relative Expression (2^-ddCt)'] - gene_data['Fold Change Min']
    gene_data['Error_Upper'] = gene_data['Fold Change Max'] - gene_data['Relative Expression (2^-ddCt)']
    
    # Combine upper and lower into one array for matplotlib errorbar
    errors = gene_data[['Error_Lower', 'Error_Upper']].values.T

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create the bar plot with error bars
    ax.bar(
        gene_data['Sample Name'], 
        gene_data['Relative Expression (2^-ddCt)'], 
        yerr=errors, 
        capsize=5, 
        color='#1f77b4', 
        alpha=0.7
    )
    
    # Draw a line at y=1 (Reference/Control Expression Level)
    ax.axhline(1, color='red', linestyle='--', linewidth=1, label='Reference Level (1.0)')

    ax.set_title(f'Relative Expression of {gene}', fontsize=16, fontweight='bold')
    ax.set_ylabel('Relative Expression (Fold Change)', fontsize=12)
    ax.set_xlabel('Sample Group', fontsize=12)
    
    # Style improvements
    ax.grid(axis='y', linestyle=':', alpha=0.6)
    ax.set_axisbelow(True)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    return fig

# ------------------------------------------------------------
# MAIN STREAMLIT APP
# ------------------------------------------------------------

def main():
    st.title("🧬 Universal qPCR ΔΔCt Analysis Tool")
    st.markdown("Upload your raw qPCR results file (Excel or CSV) for 2^-ΔΔCt calculation and visualization.")
    
    uploaded_file = st.file_uploader("Choose a file (Excel or CSV)", type=['xlsx', 'xls', 'csv'])

    if uploaded_file:
        file_extension = uploaded_file.name.split('.')[-1].lower()
        
        # --- 1. Data Loading and Column Identification ---
        with st.spinner("Parsing file and identifying columns..."):
            all_data, col_names = parse_qpcr_data(uploaded_file, file_extension)

        if all_data is not None:
            
            # --- 2. Data Cleaning and Standardization (New Step) ---
            st.success(f"Data parsed successfully! Identified columns: Sample Name ('{col_names['sample']}'), Target Name ('{col_names['target']}'), Ct ('{col_names['ct']}')")
            
            # Standardization to handle '1' and '1.0' as the same sample group
            all_data = standardize_sample_names(all_data, col_names['sample'])
            
            # Clean and filter data
            cleaned_data = clean_data(all_data, col_names)
            
            if cleaned_data.empty:
                st.error("No valid Ct data remaining after cleaning (e.g., all Ct values were non-numeric or above 35).")
                return

            # --- 3. User Selection and Parameters ---
            
            available_targets = sorted(cleaned_data['Target Name'].unique().tolist())
            available_samples = sorted(cleaned_data['Sample Name'].unique().tolist())

            st.sidebar.header("Analysis Parameters")

            # Selection for Housekeeping Gene
            hk_gene = st.sidebar.selectbox(
                "1. Select Housekeeping Gene (Reference Gene)",
                options=[g for g in available_targets if g != 'NTC'],
                index=0 if len(available_targets) > 0 and available_targets[0] != 'NTC' else 0,
                help="The gene used for normalization (ΔCt calculation)."
            )

            # Selection for Reference Group
            ref_group = st.sidebar.selectbox(
                "2. Select Reference Sample Group (Control)",
                options=available_samples,
                index=0,
                help="The sample group used as the baseline for fold change (ΔΔCt calculation). Fold change will be 1.0 for this group."
            )
            
            # List of genes to analyze (excluding HK gene)
            target_genes = [g for g in available_targets if g != hk_gene and g != 'NTC']
            
            # Selection for genes to plot
            st.sidebar.markdown("---")
            selected_genes = st.sidebar.multiselect(
                "3. Select Genes to Plot (Excluding HK)",
                options=target_genes,
                default=target_genes,
                help="Select which target genes (excluding your housekeeping gene) you want to visualize."
            )

            # --- 4. Calculation and Results Display ---
            if st.sidebar.button("Run qPCR Analysis"):
                
                # Filter data to only include the HK gene and the genes selected for analysis
                targets_to_include = selected_genes + [hk_gene]
                analysis_data = cleaned_data[cleaned_data['Target Name'].isin(targets_to_include)]
                
                if analysis_data.empty:
                    st.warning("Selected genes or housekeeping gene not found in the cleaned data.")
                    return
                
                with st.spinner("Performing 2^-ΔΔCt calculation..."):
                    summary = calculate_ddct(analysis_data, hk_gene, ref_group)
                    
                    st.subheader("Results: Relative Gene Expression (2^-ΔΔCt)")
                    
                    # Formatting the table for display
                    summary_display = summary.copy()
                    summary_display = summary_display.round({
                        'dCt_Sample_Mean': 3,
                        'dCt_SE': 3,
                        'ddCt': 3,
                        'Relative Expression (2^-ddCt)': 3,
                        'Fold Change Min': 3,
                        'Fold Change Max': 3
                    })
                    
                    # Style the reference group to have a relative expression of 1.0 (since it is the reference)
                    def highlight_ref(s):
                        if s['Sample Name'] == ref_group:
                            return ['background-color: #e6ffe6'] * len(s) # Light green background
                        return [''] * len(s)

                    st.dataframe(
                        summary_display.style.apply(highlight_ref, axis=1),
                        use_container_width=True
                    )

                    # Download button for the results
                    csv_buffer = summary.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="⬇️ Download Full Results (CSV)",
                        data=csv_buffer,
                        file_name=f"qpcr_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                        mime="text/csv"
                    )
                    
                    st.markdown("---")
                    
                    # --- 5. Graph Generation ---
                    st.subheader("Relative Expression Graphs")

                    if selected_genes:
                        with st.spinner("Creating graphs..."):
                            
                            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
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
                        st.warning("Please select at least one gene to generate graphs.")
                
            st.markdown("---")
            st.subheader("Raw Data Preview")
            st.dataframe(cleaned_data, use_container_width=True)


if __name__ == "__main__":
    main()
