"""
qpcr_analysis_app.py — Complete qPCR Analysis Web Application
------------------------------------------------------------
Features:
  - User-driven file upload for Excel (.xlsx) or CSV (.csv) files.
  - Handles multiple files for combined analysis.
  - Flexible parsing: scans data frames for Sample Name, Target Name, Ct columns.
  - User selects data start row, control gene, and reference group.
  - Calculates relative gene expression (2^-ΔΔCt) with robust statistics.
  - Generates downloadable results and bar graphs.

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
    page_title="Dynamic qPCR Analysis Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ------------------------------------------------------------
# FLEXIBLE DATA SCANNING FUNCTIONS
# ------------------------------------------------------------

@st.cache_data
def scan_for_columns(df):
    """
    Scan DataFrame for Sample Name, Target Name, and Ct columns.
    Returns dict with column names if found, None otherwise.
    """
    
    patterns = {
        'sample': ['sample name', 'sample', 'samplename', 'name'],
        'target': ['target name', 'target', 'targetname', 'gene', 'assay', 'reporter'],
        'ct': [
            'ct', 'c t', 'c-t', 'cq', 'value', 'quantification cycle',
            'cт', 'c т', 'ст', # Cyrillic variations
            'ct/cт/cq/value', 'ct value', 'c t value', 'cq value',
            'ct mean', 'cт mean',
        ]
    }
    
    def normalize_header(header):
        return str(header).strip().lower().replace(' ', '').replace('-', '').replace('/', '')

    df_cols = [normalize_header(col) for col in df.columns]
    found_cols = {}
    
    for col_type, variations in patterns.items():
        normalized_variations = [normalize_header(v) for v in variations]

        # Find the first column header that matches any of the defined normalized variations
        found_index = next(
            (i for i, df_col in enumerate(df_cols) if df_col in normalized_variations),
            None
        )
        
        if found_index is not None:
            found_cols[col_type] = df.columns[found_index]
        else:
            return None 

    return found_cols

@st.cache_data(show_spinner=False)
def load_and_preprocess_data(uploaded_files, header_rows):
    """
    Loads data from uploaded file objects using the user-specified header row.
    Combines, cleans, and prepares data for DDCT calculation.
    """
    st.info("Starting data parsing and combination...")
    all_data = []

    for i, uploaded_file in enumerate(uploaded_files):
        file_name = uploaded_file.name
        # header_rows is 1-based, pandas 'header' is 0-based index
        header_row_index = header_rows[i] - 1
        
        try:
            # Determine file type and load
            if file_name.endswith('.xlsx'):
                # For Excel, we read only the first sheet (or specify sheet_name=0)
                df = pd.read_excel(uploaded_file, header=header_row_index, engine='openpyxl')
            elif file_name.endswith('.csv'):
                df = pd.read_csv(uploaded_file, header=header_row_index)
            else:
                st.warning(f"Skipping file **{file_name}**: Unsupported format. Please upload .xlsx or .csv.")
                continue

            # Strip spaces from column names to help matching logic
            df.columns = [col.strip() for col in df.columns]

            col_map = scan_for_columns(df)
            
            if col_map:
                df_clean = df.rename(columns={
                    col_map['sample']: 'Sample Name',
                    col_map['target']: 'Target Name',
                    col_map['ct']: 'Ct'
                })
                # Ensure Ct is numeric, coercing non-numeric (e.g., 'Undetermined') to NaN
                df_clean['Ct'] = pd.to_numeric(df_clean['Ct'], errors='coerce')
                # Keep only essential columns and drop rows where Ct is NaN
                df_clean = df_clean[['Sample Name', 'Target Name', 'Ct']].dropna(subset=['Ct'])
                
                if not df_clean.empty:
                    all_data.append(df_clean)
                    st.success(f"✅ Processed file: **{file_name}**")
                else:
                    st.warning(f"⚠️ File **{file_name}** had no valid Ct data after cleaning. Skipped.")
            else:
                st.error(f"❌ Failed to find mandatory columns (Sample Name, Target Name, Ct) in file **{file_name}**. Check your header row setting.")

        except Exception as e:
            st.error(f"Error loading file **{file_name}**: {e}")

    if not all_data:
        st.error("No valid qPCR data could be loaded for analysis.")
        return None, None, None, None

    # Combine all valid dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Calculate group name for each sample (assuming Sample Name is the unique group identifier)
    combined_df['Group'] = combined_df['Sample Name'].apply(lambda x: str(x).strip())
    
    samples = combined_df['Sample Name'].unique().tolist()
    targets = combined_df['Target Name'].unique().tolist()
    groups = combined_df['Group'].unique().tolist()
    
    st.success(f"Data loading complete. Found {len(samples)} unique samples and {len(targets)} unique targets.")
    
    return combined_df, samples, targets, groups

# ------------------------------------------------------------
# CALCULATION FUNCTIONS (Unchanged from previous versions logic)
# ------------------------------------------------------------

def calculate_ddct(data, control_gene, reference_group):
    """Performs the 2^-ΔΔCt calculation and T-Test for significance."""
    
    # 1. Group by Sample and Target, calculate Mean Ct
    ct_mean = data.groupby(['Group', 'Target Name'])['Ct'].mean().reset_index()
    ct_mean = ct_mean.rename(columns={'Ct': 'Mean Ct'})
    
    # 2. Calculate ΔCt: (Target Mean Ct) - (Control Gene Mean Ct)
    control_cts = ct_mean[ct_mean['Target Name'] == control_gene].set_index('Group')['Mean Ct']
    if control_cts.empty:
        st.error(f"The selected control gene '{control_gene}' was not found in the data.")
        return None
    
    ct_mean['Control Ct'] = ct_mean['Group'].map(control_cts)
    ct_mean = ct_mean.dropna(subset=['Control Ct']).copy() 
    ct_mean['ΔCt'] = ct_mean['Mean Ct'] - ct_mean['Control Ct']

    # 3. Calculate ΔΔCt: (ΔCt of Sample) - (ΔCt of Reference Group)
    ref_dct = ct_mean[
        (ct_mean['Group'] == reference_group) & 
        (ct_mean['Target Name'] != control_gene)
    ].set_index('Target Name')['ΔCt']
    
    if ref_dct.empty:
        st.error(f"The reference group '{reference_group}' did not have detectable ΔCt values.")
        return None
        
    def get_ref_dct(row):
        return ref_dct.get(row['Target Name'])

    ct_mean['Reference ΔCt'] = ct_mean.apply(get_ref_dct, axis=1)
    summary_df = ct_mean[ct_mean['Target Name'] != control_gene].copy()
    summary_df = summary_df.dropna(subset=['Reference ΔCt']).copy()
    summary_df['ΔΔCt'] = summary_df['ΔCt'] - summary_df['Reference ΔCt']
    
    # 4. Calculate Relative Expression (Fold Change): 2^-ΔΔCt
    summary_df['Fold Change (2^-ΔΔCt)'] = 2**(-summary_df['ΔΔCt'])
    
    # 5. Calculate Standard Error for ΔΔCt
    raw_data_agg = data.groupby(['Group', 'Target Name'])['Ct'].agg(['mean', 'std', 'count']).reset_index()
    raw_data_agg = raw_data_agg.rename(columns={'mean': 'Mean Ct', 'std': 'Std Ct', 'count': 'N'})
    
    control_agg = raw_data_agg[raw_data_agg['Target Name'] == control_gene].set_index('Group')[['Std Ct', 'N']]
    control_agg = control_agg.rename(columns={'Std Ct': 'Std Ct Control', 'N': 'N Control'})
    
    full_agg = raw_data_agg.merge(control_agg, on='Group', how='left')
    full_agg = full_agg[full_agg['Target Name'] != control_gene].dropna(subset=['Std Ct Control']).copy()

    # SE(ΔCt) calculation
    full_agg['SE ΔCt'] = np.sqrt(
        (full_agg['Std Ct']**2 / full_agg['N']) + (full_agg['Std Ct Control']**2 / full_agg['N Control'])
    )
    
    ref_se_dct = full_agg[full_agg['Group'] == reference_group].set_index('Target Name')['SE ΔCt']
    full_agg['Ref SE ΔCt'] = full_agg['Target Name'].map(ref_se_dct)
    full_agg = full_agg.dropna(subset=['Ref SE ΔCt']).copy()
    
    # SE ΔΔCt calculation
    full_agg['SE ΔΔCt'] = np.sqrt(full_agg['SE ΔCt']**2 + full_agg['Ref SE ΔCt']**2)
    
    summary_df_se = summary_df.merge(full_agg[['Group', 'Target Name', 'SE ΔΔCt']], 
                                      on=['Group', 'Target Name'], how='left')
    
    summary_df_se = summary_df_se.set_index(['Group', 'Target Name'])
    
    # Standard Deviation of Fold Change (Upper/Lower Bound for Error Bars)
    summary_df_se['Upper Bound'] = 2**(-(summary_df_se['ΔΔCt'] - summary_df_se['SE ΔΔCt']))
    summary_df_se['Lower Bound'] = 2**(-(summary_df_se['ΔΔCt'] + summary_df_se['SE ΔΔCt']))
    
    # 6. T-Test for Statistical Significance (comparing ΔCt of Sample vs ΔCt of Reference)
    p_values = {}
    for (group, target), _ in summary_df_se.iterrows():
        if group == reference_group:
            p_values[(group, target)] = np.nan
            continue
            
        group_dct_raw = []
        for sample in data[data['Group'] == group]['Sample Name'].unique():
            target_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == target)]['Ct'].mean()
            control_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == control_gene)]['Ct'].mean()
            if not pd.isna(target_ct) and not pd.isna(control_ct):
                group_dct_raw.append(target_ct - control_ct)
                
        ref_dct_raw = []
        for sample in data[data['Group'] == reference_group]['Sample Name'].unique():
            target_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == target)]['Ct'].mean()
            control_ct = data[(data['Sample Name'] == sample) & (data['Target Name'] == control_gene)]['Ct'].mean()
            if not pd.isna(target_ct) and not pd.isna(control_ct):
                ref_dct_raw.append(target_ct - control_ct)

        if len(group_dct_raw) > 1 and len(ref_dct_raw) > 1:
            t_stat, p_val = stats.ttest_ind(group_dct_raw, ref_dct_raw, equal_var=False)
            p_values[(group, target)] = p_val
        else:
            p_values[(group, target)] = np.nan
            
    summary_df_se['P-Value'] = pd.Series(p_values)
    
    # 7. Final Formatting and Cleanup
    summary_df_se = summary_df_se.reset_index()
    summary_df_se['Significant'] = summary_df_se['P-Value'].apply(
        lambda x: 'Yes (p < 0.05)' if x < 0.05 else ('No' if not pd.isna(x) else 'N/A')
    )
    
    final_cols = [
        'Group', 'Target Name', 'Mean Ct', 'ΔCt', 'ΔΔCt', 
        'Fold Change (2^-ΔΔCt)', 'SE ΔΔCt', 'Upper Bound', 'Lower Bound', 'P-Value', 'Significant'
    ]
    
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
        
    ref_group = plot_data[plot_data['Fold Change (2^-ΔΔCt)'] == 1]['Group'].iloc[0] if not plot_data[plot_data['Fold Change (2^-ΔΔCt)'] == 1].empty else None
    
    # Calculate error bar height (distance from Fold Change to Upper Bound)
    plot_data['Error'] = plot_data['Upper Bound'] - plot_data['Fold Change (2^-ΔΔCt)']
    plot_data.loc[plot_data['Group'] == ref_group, 'Error'] = 0.0
    
    plot_data['Color'] = plot_data['Significant'].apply(
        lambda x: 'darkred' if 'Yes' in x else ('cornflowerblue' if 'No' in x else 'gray')
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.bar(
        plot_data['Group'], 
        plot_data['Fold Change (2^-ΔΔCt)'], 
        yerr=plot_data['Error'], 
        capsize=5, 
        color=plot_data['Color']
    )

    # Add significance markers
    for i, (idx, row) in enumerate(plot_data.iterrows()):
        if 'Yes' in row['Significant']:
            ax.text(i, row['Upper Bound'] + 0.1, '*', ha='center', va='bottom', fontsize=18, color='darkred')
        
    ax.axhline(1, color='gray', linestyle='--', label='Reference (Fold Change = 1)')
    ax.set_ylabel('Fold Change (2⁻ΔΔCt) Relative to Reference')
    ax.set_xlabel('Sample Group')
    ax.set_title(f'Relative Gene Expression: {gene}')
    ax.set_ylim(0, plot_data['Upper Bound'].max() * 1.25)
    
    # Improve aesthetics
    plt.xticks(rotation=45, ha='right')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    
    return fig

# ------------------------------------------------------------
# STREAMLIT MAIN APP
# ------------------------------------------------------------

def main():
    """Main function for the Streamlit app layout and logic."""
    
    st.title("🧬 Dynamic qPCR 2⁻ΔΔCt Analysis Tool")
    st.markdown("Upload your raw qPCR result files (.xlsx or .csv) to begin the analysis.")
    
    st.sidebar.header("1. Upload Data")
    
    # 1. File Uploader
    uploaded_files = st.sidebar.file_uploader(
        "Upload one or more qPCR results files",
        type=['csv', 'xlsx'],
        accept_multiple_files=True
    )
    
    if not uploaded_files:
        st.info("Please upload your raw data files in the sidebar to proceed.")
        return

    # 2. Header Row Selector
    st.sidebar.markdown("---")
    st.sidebar.header("2. Configure Data Files")
    
    header_rows = []
    
    # Ensure a Header Row is selected for each file, defaulting to row 1 (for files with no metadata)
    # or allowing the user to set the correct row (e.g., row 37 for some exports).
    for i, file in enumerate(uploaded_files):
        header_row = st.sidebar.number_input(
            f"Header Row (1-based index) for: **{file.name}**",
            min_value=1, 
            value=1, # Default to 1 (A1 header)
            step=1, 
            key=f'header_{i}'
        )
        header_rows.append(header_row)

    # 3. Load and Preprocess Data
    data, samples, targets, groups = load_and_preprocess_data(uploaded_files, header_rows)
    
    if data is None:
        # Stop execution if data loading failed
        return 
        
    st.sidebar.header("3. Analysis Parameters")
    
    # 4. Select Housekeeping/Control Gene
    default_control_index = targets.index('ACTIN') if 'ACTIN' in targets else (0 if targets else None)
    control_gene = st.sidebar.selectbox(
        "Select Housekeeping/Control Gene:",
        targets,
        index=default_control_index if default_control_index is not None else 0
    )
    
    # 5. Select Reference Group
    reference_group = st.sidebar.selectbox(
        "Select Reference Group (Fold Change = 1):",
        groups,
        index=0
    )
    
    # 6. Calculate Button
    if st.sidebar.button("🔬 Start 2⁻ΔΔCt Calculation"):
        summary = calculate_ddct(data, control_gene, reference_group)
        
        if summary is not None and not summary.empty:
            st.session_state['summary'] = summary
            st.session_state['targets'] = targets
            st.session_state['control_gene'] = control_gene
            st.success("Calculation complete! Scroll down for results and visualizations.")
        
    # --- Display Results ---
    if 'summary' in st.session_state:
        summary = st.session_state['summary']
        targets = st.session_state['targets']
        control_gene = st.session_state['control_gene']

        st.markdown("---")
        st.subheader("4. Results: Relative Gene Expression Summary")
        
        # Display summary table
        display_summary = summary.copy()
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
        st.subheader("5. Visualization & Graphs")
        
        target_genes = [t for t in targets if t != control_gene]
        selected_genes = st.multiselect(
            "Select target genes to plot (Fold Change vs. Group):",
            target_genes,
            default=target_genes
        )
        
        if st.button("📈 Generate Graphs"):
            if selected_genes:
                with st.spinner("Creating graphs..."):
                    
                    zip_buffer = io.BytesIO()
                    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                        
                        for gene in selected_genes:
                            fig = create_bar_graph(summary, gene)
                            if fig:
                                # Save the figure to a buffer as PNG
                                img_buffer = io.BytesIO()
                                fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                                img_buffer.seek(0)
                                
                                # Add the PNG to the zip file
                                zip_file.writestr(f"{gene}_expression.png", img_buffer.getvalue())
                                
                                # Display the figure in Streamlit
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


if __name__ == "__main__":
    main()
