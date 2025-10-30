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
    # NOTE: All DataFrame headers are lowercased and stripped of leading/trailing space 
    # before comparison against this list, to ensure flexible matching.
    patterns = {
        'sample': ['sample name', 'sample', 'samplename', 'name'],
        'target': ['target name', 'target', 'targetname', 'gene', 'assay', 'reporter'],
        'ct': [
            # Standard variations
            'ct', 'c t', 'c-t', 'cq', 'value', 'quantification cycle',

            # Cyrillic/Mixed variations (e.g., from StepOne software exports)
            'cт',  # Lowercase Cyrillic 'te' (Cт)
            'c т', # Cyrillic with space
            'ст',  # Lowercase Cyrillic 'es-te' (less common, but included for robustness)
            
            # Variations with slashes, words, and spaces often seen in raw files
            'ct/cт/cq/value', 
            'ct value', 'c t value', 'cq value',
            'ct mean', 'cт mean',
        ]
    }
    
    df_cols = [col.strip().lower() for col in df.columns]
    
    # Mapping found column type (e.g., 'sample') to its actual column name in the DataFrame
    found_cols = {}
    
    for col_type, variations in patterns.items():
        # Find the first column header that matches any of the defined variations
        found_name = next(
            (df.columns[i] for i, df_col in enumerate(df_cols) if df_col in variations),
            None
        )
        
        if found_name:
            found_cols[col_type] = found_name
        else:
            # If any mandatory column is missing, the scan fails for this sheet
            return None 

    # Check for empty columns (e.g., if a column header was found but data is missing)
    for col_name in found_cols.values():
        if df[col_name].isnull().all() or not pd.api.types.is_numeric_dtype(pd.to_numeric(df[col_name], errors='coerce').dropna()):
            # If a CT column is found but non-numeric, it's not the correct data column
            if col_name == found_cols.get('ct'):
                 return None
            
    return found_cols

@st.cache_data
def load_and_preprocess_data(uploaded_file):
    """Loads uploaded file, scans all sheets/data frames, and combines results."""
    
    st.info("Scanning file for Sample, Target, and Ct columns across all sheets/data...")
    all_data = []

    try:
        # Try loading as Excel (xlsx, xls)
        xls = pd.ExcelFile(uploaded_file, engine='openpyxl')
        sheet_names = xls.sheetnames
        
        # Iterate over all sheets
        for sheet_name in sheet_names:
            df = pd.read_excel(xls, sheet_name=sheet_name)
            col_map = scan_for_columns(df)
            
            if col_map:
                df_clean = df.rename(columns={
                    col_map['sample']: 'Sample Name',
                    col_map['target']: 'Target Name',
                    col_map['ct']: 'Ct'
                })
                # Ensure Ct is numeric, coercing non-numeric (e.g., 'Undetermined') to NaN
                df_clean['Ct'] = pd.to_numeric(df_clean['Ct'], errors='coerce')
                df_clean = df_clean[['Sample Name', 'Target Name', 'Ct']].dropna(subset=['Ct'])
                
                if not df_clean.empty:
                    all_data.append(df_clean)
                    st.success(f"Successfully processed sheet/data: **{sheet_name}**")
            else:
                st.warning(f"Skipped sheet/data: **{sheet_name}** (Required columns not found or Ct column is empty/non-numeric)")

    except Exception as e_excel:
        # If Excel load failed, try CSV (or other single-sheet formats)
        try:
            uploaded_file.seek(0)
            df = pd.read_csv(uploaded_file)
            col_map = scan_for_columns(df)
            
            if col_map:
                df_clean = df.rename(columns={
                    col_map['sample']: 'Sample Name',
                    col_map['target']: 'Target Name',
                    col_map['ct']: 'Ct'
                })
                df_clean['Ct'] = pd.to_numeric(df_clean['Ct'], errors='coerce')
                df_clean = df_clean[['Sample Name', 'Target Name', 'Ct']].dropna(subset=['Ct'])
                
                if not df_clean.empty:
                    all_data.append(df_clean)
                    st.success("Successfully processed CSV/single-data file.")
                else:
                    st.error("Data found but became empty after filtering out non-numeric Ct values.")
            else:
                st.error("Required columns not found in the file.")
                
        except Exception as e_csv:
            st.error(f"Error loading file: Could not read as Excel ({e_excel}) or CSV ({e_csv}).")
            return None, None, None, None

    if not all_data:
        st.error("No valid qPCR data found in the uploaded file. Please check your data format.")
        return None, None, None, None

    # Combine all valid dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Calculate group name for each sample
    combined_df['Group'] = combined_df['Sample Name'].apply(lambda x: str(x).strip())
    
    samples = combined_df['Sample Name'].unique().tolist()
    targets = combined_df['Target Name'].unique().tolist()
    groups = combined_df['Group'].unique().tolist()
    
    st.success(f"Data loading complete. Found {len(samples)} unique samples and {len(targets)} unique targets.")
    
    return combined_df, samples, targets, groups

# ------------------------------------------------------------
# CALCULATION FUNCTIONS
# ------------------------------------------------------------

def calculate_ddct(data, control_gene, reference_group):
    """Performs the 2^-ΔΔCt calculation."""
    
    # 1. Group by Sample and Target, calculate Mean Ct
    ct_mean = data.groupby(['Group', 'Target Name'])['Ct'].mean().reset_index()
    ct_mean = ct_mean.rename(columns={'Ct': 'Mean Ct'})
    
    # 2. Calculate ΔCt: (Target Mean Ct) - (Control Gene Mean Ct)
    control_cts = ct_mean[ct_mean['Target Name'] == control_gene].set_index('Group')['Mean Ct']
    
    if control_cts.empty:
        st.error(f"The selected control gene '{control_gene}' was not found in the data.")
        return None
    
    # Map control Ct back to the main dataframe by 'Group'
    ct_mean['Control Ct'] = ct_mean['Group'].map(control_cts)
    
    # Filter out samples where the control gene was not detected
    ct_mean = ct_mean.dropna(subset=['Control Ct']).copy() 
    
    ct_mean['ΔCt'] = ct_mean['Mean Ct'] - ct_mean['Control Ct']

    # 3. Calculate ΔΔCt: (ΔCt of Sample) - (ΔCt of Reference Group)
    ref_dct = ct_mean[
        (ct_mean['Group'] == reference_group) & 
        (ct_mean['Target Name'] != control_gene) # Don't use control gene itself for ref_dct
    ].set_index('Target Name')['ΔCt']
    
    if ref_dct.empty:
        st.error(f"The reference group '{reference_group}' did not have detectable ΔCt values.")
        return None
        
    # Map reference ΔCt back to the main dataframe by 'Target Name'
    def get_ref_dct(row):
        return ref_dct.get(row['Target Name'])

    ct_mean['Reference ΔCt'] = ct_mean.apply(get_ref_dct, axis=1)

    # Filter out the control gene rows, as they don't have a meaningful ΔΔCt (they would be 0 or NaN)
    summary_df = ct_mean[ct_mean['Target Name'] != control_gene].copy()
    
    # Filter out rows where the reference ΔCt could not be found
    summary_df = summary_df.dropna(subset=['Reference ΔCt']).copy()

    summary_df['ΔΔCt'] = summary_df['ΔCt'] - summary_df['Reference ΔCt']
    
    # 4. Calculate Relative Expression (Fold Change): 2^-ΔΔCt
    summary_df['Fold Change (2^-ΔΔCt)'] = 2**(-summary_df['ΔΔCt'])
    
    # 5. Calculate Standard Error for ΔΔCt
    # This requires looking back at the raw Ct data for replicate variance.
    summary_df = summary_df.set_index(['Group', 'Target Name'])
    
    # Get raw Ct replicates
    raw_ct_data = data.merge(summary_df.reset_index()[['Group', 'Target Name']], 
                             on=['Group', 'Target Name'], how='inner')
                             
    # Calculate ΔCt for each raw replicate (Ct_target - Ct_control)
    def calculate_raw_dct(group):
        target_cts = group[group['Target Name'] == group['Target Name'].iloc[0]]['Ct'].values
        control_cts = raw_ct_data[(raw_ct_data['Group'] == group['Group'].iloc[0]) & 
                                  (raw_ct_data['Target Name'] == control_gene)]['Ct'].values
                                  
        if len(control_cts) == 0:
            return np.nan
            
        # Perform all possible replicate-to-replicate subtractions (only if replicate numbers match)
        # Assuming triplicates (or equal number of replicates for target and control)
        # If not matched, we use the average of control_cts for safety, but this loses variance info.
        # For simplicity and robustness in diverse raw data, we stick to the primary summary df:
        
        # Get the standard deviation of the raw ΔCt replicates for the current group/target
        raw_target = data[(data['Group'] == group['Group'].iloc[0]) & (data['Target Name'] == group['Target Name'].iloc[0])]['Ct']
        raw_control = data[(data['Group'] == group['Group'].iloc[0]) & (data['Target Name'] == control_gene)]['Ct']
        
        # Simple approach for SE calculation (variance addition, then SE of mean):
        # Var(ΔCt) = Var(Ct_Target) + Var(Ct_Control)
        var_dct = raw_target.var() + raw_control.var()

        # StdDev(ΔCt) = sqrt(Var(ΔCt))
        # N_Target and N_Control are assumed to be equal for the same group (e.g., n=3)
        n_replicates = raw_target.count()
        
        # We need SE for the final ΔΔCt, which is sqrt(SE(ΔCt_Sample)^2 + SE(ΔCt_Ref)^2)
        
        # Let's calculate SE of ΔCt for all groups first
        raw_data_agg = data.groupby(['Group', 'Target Name'])['Ct'].agg(['mean', 'std', 'count']).reset_index()
        raw_data_agg = raw_data_agg.rename(columns={'mean': 'Mean Ct', 'std': 'Std Ct', 'count': 'N'})
        
        # Add control gene Std Ct
        control_agg = raw_data_agg[raw_data_agg['Target Name'] == control_gene].set_index('Group')[['Std Ct', 'N']]
        control_agg = control_agg.rename(columns={'Std Ct': 'Std Ct Control', 'N': 'N Control'})
        
        full_agg = raw_data_agg.merge(control_agg, on='Group', how='left')
        
        # Filter out control gene rows and ensure control data is present
        full_agg = full_agg[full_agg['Target Name'] != control_gene].dropna(subset=['Std Ct Control']).copy()

        # Variance of ΔCt = Var(Target) + Var(Control) = Std^2 + Std Control^2
        full_agg['Var ΔCt'] = full_agg['Std Ct']**2 + full_agg['Std Ct Control']**2
        
        # Standard Error of ΔCt = sqrt(Var(ΔCt))
        # No, SE(ΔCt) = sqrt( (Std(Ct_T)^2 / N_T) + (Std(Ct_C)^2 / N_C) )
        full_agg['SE ΔCt'] = np.sqrt(
            (full_agg['Std Ct']**2 / full_agg['N']) + (full_agg['Std Ct Control']**2 / full_agg['N Control'])
        )
        
        # Reference SE ΔCt
        ref_se_dct = full_agg[full_agg['Group'] == reference_group].set_index('Target Name')['SE ΔCt']
        full_agg['Ref SE ΔCt'] = full_agg['Target Name'].map(ref_se_dct)
        
        # Filter out rows where reference SE is missing
        full_agg = full_agg.dropna(subset=['Ref SE ΔCt']).copy()
        
        # SE ΔΔCt = sqrt(SE(ΔCt_Sample)^2 + SE(ΔCt_Ref)^2)
        full_agg['SE ΔΔCt'] = np.sqrt(full_agg['SE ΔCt']**2 + full_agg['Ref SE ΔCt']**2)
        
        # Merge SE ΔΔCt back into the summary_df
        summary_df_se = summary_df.merge(full_agg[['Group', 'Target Name', 'SE ΔΔCt']], 
                                          on=['Group', 'Target Name'], how='left')
        
        summary_df_se = summary_df_se.set_index(['Group', 'Target Name'])
        
        # Standard Deviation of Fold Change (2^-ΔΔCt)
        summary_df_se['Upper Bound'] = 2**(-(summary_df_se['ΔΔCt'] - summary_df_se['SE ΔΔCt']))
        summary_df_se['Lower Bound'] = 2**(-(summary_df_se['ΔΔCt'] + summary_df_se['SE ΔΔCt']))
        
        # 6. T-Test for Statistical Significance (comparing ΔCt of Sample vs ΔCt of Reference)
        p_values = {}
        for (group, target), row in summary_df_se.iterrows():
            if group == reference_group:
                p_values[(group, target)] = np.nan # Reference group is always 1, not tested
                continue
                
            # Get raw ΔCt values for this sample group
            group_dct_raw = []
            for sample in data[data['Group'] == group]['Sample Name'].unique():
                target_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == target)]['Ct'].mean()
                control_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == control_gene)]['Ct'].mean()
                if not pd.isna(target_ct) and not pd.isna(control_ct):
                    group_dct_raw.append(target_ct - control_ct)
                    
            # Get raw ΔCt values for reference group
            ref_dct_raw = []
            for sample in data[data['Group'] == reference_group]['Sample Name'].unique():
                target_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == target)]['Ct'].mean()
                control_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == control_gene)]['Ct'].mean()
                if not pd.isna(target_ct) and not pd.isna(control_ct):
                    ref_dct_raw.append(target_ct - control_ct)

            if len(group_dct_raw) > 1 and len(ref_dct_raw) > 1:
                # Perform two-sided T-Test assuming unequal variance (Welch's T-Test)
                t_stat, p_val = stats.ttest_ind(group_dct_raw, ref_dct_raw, equal_var=False)
                p_values[(group, target)] = p_val
            else:
                p_values[(group, target)] = np.nan
                
        summary_df_se['P-Value'] = pd.Series(p_values)
        
        # 7. Final Formatting and Cleanup
        summary_df_se = summary_df_se.reset_index()
        summary_df_se['Significant'] = summary_df_se['P-Value'].apply(
            lambda x: 'Yes' if x < 0.05 else ('No' if not pd.isna(x) else 'N/A')
        )
        
        final_cols = [
            'Group', 'Target Name', 'Mean Ct', 'ΔCt', 'ΔΔCt', 
            'Fold Change (2^-ΔΔCt)', 'SE ΔΔCt', 'Upper Bound', 'Lower Bound', 'P-Value', 'Significant'
        ]
        
        # Drop columns that were used for intermediate calculation but might not exist if previous steps failed
        final_summary = summary_df_se.filter(items=final_cols, axis=1)

        return final_summary

# ------------------------------------------------------------
# PLOTTING FUNCTIONS
# ------------------------------------------------------------

def create_bar_graph(summary, gene):
    """Creates a bar plot for the specified gene across all groups."""
    
    plot_data = summary[summary['Target Name'] == gene].copy()
    
    if plot_data.empty:
        return None
        
    # Exclude reference group from SE error bar calculation for clarity (SE at 1 is usually 0)
    ref_group = plot_data[plot_data['Fold Change (2^-ΔΔCt)'] == 1]['Group'].iloc[0] if not plot_data[plot_data['Fold Change (2^-ΔΔCt)'] == 1].empty else None
    
    # Calculate error bar height (distance from Fold Change to Upper Bound)
    plot_data['Error'] = plot_data['Upper Bound'] - plot_data['Fold Change (2^-ΔΔCt)']
    plot_data.loc[plot_data['Group'] == ref_group, 'Error'] = 0.0 # Ref group error is 0
    
    # Simple check for significance for color coding
    plot_data['Color'] = plot_data['Significant'].apply(lambda x: 'darkred' if x == 'Yes' else ('cornflowerblue' if x == 'No' else 'gray'))

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Use the calculated error for error bars
    bars = ax.bar(
        plot_data['Group'], 
        plot_data['Fold Change (2^-ΔΔCt)'], 
        yerr=plot_data['Error'], 
        capsize=5, 
        color=plot_data['Color']
    )

    # Label significance above bars
    for i, row in plot_data.iterrows():
        if row['Significant'] == 'Yes':
            ax.text(i, row['Upper Bound'] + 0.1, '*', ha='center', va='bottom', fontsize=18, color='darkred')
        
    ax.axhline(1, color='gray', linestyle='--')
    ax.set_ylabel('Fold Change (2⁻ΔΔCt) Relative to Reference')
    ax.set_xlabel('Sample Group')
    ax.set_title(f'Relative Gene Expression: {gene}')
    ax.set_ylim(0, plot_data['Upper Bound'].max() * 1.2)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    return fig

# ------------------------------------------------------------
# STREAMLIT MAIN APP
# ------------------------------------------------------------

def main():
    """Main function for the Streamlit app layout and logic."""
    
    st.title("🧬 qPCR 2⁻ΔΔCt Analysis Tool")
    
    # Load Firebase config (required for the canvas environment)
    firebaseConfig = {}
    try:
        if '__firebase_config' in globals():
            firebaseConfig = JSON.parse(__firebase_config)
    except:
        pass # Ignore if not running in the specific environment

    st.sidebar.header("1. Upload Data")
    uploaded_file = st.sidebar.file_uploader(
        "Upload your qPCR results file (.xlsx or .csv)", 
        type=['xlsx', 'csv']
    )
    
    if uploaded_file is None:
        st.info("Upload a file to begin the qPCR analysis.")
        return

    # Load and preprocess data (cached)
    data, samples, targets, groups = load_and_preprocess_data(uploaded_file)
    
    if data is None:
        return # Stop if data loading failed
        
    st.sidebar.header("2. Analysis Parameters")
    
    # Step 2: Select housekeeping/control gene
    control_gene = st.sidebar.selectbox(
        "Select Housekeeping/Control Gene:",
        targets,
        index=targets.index('ACTIN') if 'ACTIN' in targets else 0
    )
    
    # Step 3: Select reference group
    # Note: Reference group is determined from the 'Sample Name' column unique values
    reference_group = st.sidebar.selectbox(
        "Select Reference Group (Fold Change = 1):",
        groups,
        index=0
    )
    
    if st.sidebar.button("Calculate 2⁻ΔΔCt"):
        summary = calculate_ddct(data, control_gene, reference_group)
        
        if summary is not None and not summary.empty:
            st.session_state['summary'] = summary
            st.session_state['targets'] = targets
            st.session_state['control_gene'] = control_gene
            st.success("Calculation complete! Results are shown below.")
        
    if 'summary' in st.session_state:
        summary = st.session_state['summary']
        targets = st.session_state['targets']
        control_gene = st.session_state['control_gene']

        st.subheader("Results: Relative Gene Expression Summary")
        
        # Display summary table
        display_summary = summary.copy()
        
        # Format the numeric columns for display
        numeric_cols = ['Mean Ct', 'ΔCt', 'ΔΔCt', 'Fold Change (2^-ΔΔCt)', 'SE ΔΔCt', 'Upper Bound', 'Lower Bound', 'P-Value']
        for col in numeric_cols:
            if col in display_summary.columns:
                display_summary[col] = display_summary[col].round(3)
        
        st.dataframe(
            display_summary,
            hide_index=True,
            use_container_width=True
        )

        # Download button for the results table
        csv_buffer = display_summary.to_csv(index=False).encode('utf-8')
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        st.download_button(
            label="⬇️ Download Results Table (CSV)",
            data=csv_buffer,
            file_name=f"qpcr_summary_{timestamp}.csv",
            mime="text/csv"
        )
        
        # --- Graphing Section ---
        st.markdown("---")
        st.subheader("3. Visualization & Graphs")
        
        target_genes = [t for t in targets if t != control_gene]
        selected_genes = st.multiselect(
            "Select target genes to plot (Fold Change vs. Group):",
            target_genes,
            default=target_genes
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
