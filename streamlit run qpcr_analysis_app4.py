"""
qpcr_analysis_app.py — Interactive and Customizable qPCR Analysis Tool
---------------------------------------------------------------------
Features:
  - Multi-file parsing: scans all sheets and multiple uploaded Excel/CSV files.
  - Robust column matching (Sample Name, Target Name, Ct).
  - Standardizes numeric sample names (e.g., '1' and '1.0' become '1').
  - INTERACTIVE DATA EDITING: Allows renaming and reordering of sample groups.
  - Calculates relative gene expression (2^-ΔΔCt).
  - Customizable bar graphs (color, Y-axis limits) for each gene.
  - Downloadable final results and graph data as an Excel file.

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
    page_title="Customizable qPCR Analysis Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ------------------------------------------------------------
# FLEXIBLE DATA SCANNING FUNCTIONS (KEPT AS IS)
# ------------------------------------------------------------

def scan_for_columns(df):
    """
    Scan DataFrame for Sample Name, Target Name, and Ct columns.
    Returns dict with column names if found, None otherwise.
    """
    patterns = {
        'sample': ['sample name', 'sample', 'samplename', 'name', 'sample id'],
        'target': ['target name', 'target', 'targetname', 'gene', 'assay', 'reporter'],
        'ct': ['ct', 'cт', 'cq', 'quantification cycle'],
    }
    
    col_map = {}
    normalized_cols = {re.sub(r'[^a-z0-9]', '', str(col).lower()): col for col in df.columns}
    
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
        return None

def parse_qpcr_data(uploaded_file, file_extension):
    """
    Loads data, finds header, concatenates all sheets, and returns a single DataFrame 
    with standardized column names ('Sample Name', 'Target Name', 'Ct').
    """
    all_data = None
    col_names = None

    try:
        if file_extension in ['xlsx', 'xls']:
            xls = pd.ExcelFile(uploaded_file)
            data_frames = []
            
            for sheet_name in xls.sheet_names:
                df = pd.read_excel(xls, sheet_name=sheet_name, header=None)
                header_row_index = -1
                
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
            
            if data_frames:
                all_data = pd.concat(data_frames, ignore_index=True)
            
        elif file_extension == 'csv':
            df = pd.read_csv(uploaded_file, header=None)
            header_row_index = -1
            
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
            
        else:
            return None

        if all_data is not None and col_names:
            standardized_data = all_data.rename(columns={
                col_names['sample']: 'Sample Name',
                col_names['target']: 'Target Name',
                col_names['ct']: 'Ct'
            })
            standardized_data = standardized_data[['Sample Name', 'Target Name', 'Ct']].copy()
            return standardized_data
        
        return None

    except Exception as e:
        st.error(f"An error occurred during parsing {uploaded_file.name}: {e}")
        return None

# ------------------------------------------------------------
# DATA CLEANING AND STANDARDIZATION (KEPT AS IS)
# ------------------------------------------------------------

def standardize_sample_names(df):
    """
    Standardizes numeric sample names (e.g., convert 1.0 to 1).
    """
    sample_col = 'Sample Name'
    s = df[sample_col].astype(str).str.strip().copy()
    numeric_s = pd.to_numeric(s, errors='coerce')
    is_numeric = numeric_s.notna()
    
    try:
        df.loc[is_numeric, sample_col] = numeric_s[is_numeric].round(0).astype(int).astype(str)
    except Exception:
        df.loc[is_numeric, sample_col] = numeric_s[is_numeric].astype(str).str.replace(r'\.0$', '', regex=True)
    
    df[sample_col] = df[sample_col].astype(str).str.strip()
    return df

def clean_data(df):
    """
    Filters out non-numeric Ct values, ensures columns are correct types, 
    and filters out high Ct values.
    """
    df['Ct'] = pd.to_numeric(df['Ct'], errors='coerce')
    df = df.dropna(subset=['Ct'])
    df = df[df['Ct'] <= 35] 
    df['Target Name'] = df['Target Name'].astype(str).str.strip()
    df['Sample Name'] = df['Sample Name'].astype(str).str.strip()
    df = df[df['Sample Name'].str.len() > 0]
    df = df[df['Target Name'].str.len() > 0]
    return df

# ------------------------------------------------------------
# NEW CALCULATION FUNCTION (REVISED)
# ------------------------------------------------------------

@st.cache_data
def calculate_ddct(df, hk_gene, ref_group):
    """
    Performs the 2^-ΔΔCt calculation using the chosen housekeeping gene and reference group.
    Uses st.cache_data for performance, rerun only if inputs change.
    """
    
    # 1. Calculate ΔCt (Ct_Target - Ct_HK) for each replicate
    # Get Housekeeping (HK) Ct for each Sample
    hk_ct = df[df['Target Name'] == hk_gene].groupby('Sample Name')['Ct'].mean().reset_index()
    hk_ct = hk_ct.rename(columns={'Ct': 'Ct_HK_Mean'})
    
    df_ddct = pd.merge(df, hk_ct, on='Sample Name', how='left')

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
    
    summary = pd.merge(dCt_summary, ref_dCt, on='Target Name', how='left')
    
    # Calculate ΔΔCt
    summary['ddCt'] = summary['dCt_Mean'] - summary['dCt_Ref_Mean']
    
    # 3. Calculate Relative Expression (2^-ΔΔCt)
    summary['Relative Expression (2^-ddCt)'] = 2**(-summary['ddCt'])
    
    # Calculate Fold Change Min and Max based on SE
    summary['Fold Change Min'] = 2**(-(summary['ddCt'] + summary['dCt_SE']))
    summary['Fold Change Max'] = 2**(-(summary['ddCt'] - summary['dCt_SE']))
    
    final_summary = summary.drop(columns=['dCt_Ref_Mean']).rename(columns={'dCt_Mean': 'dCt_Sample_Mean'})
    
    # For the reference group, ensure fold change is exactly 1.0
    final_summary.loc[final_summary['Sample Name'] == ref_group, 'Relative Expression (2^-ddCt)'] = 1.0
    
    return final_summary

# ------------------------------------------------------------
# NEW GRAPHING FUNCTION (CUSTOMIZABLE)
# ------------------------------------------------------------

def create_custom_bar_graph(summary_df, gene, sample_order, plot_options):
    """
    Creates a bar graph for a single gene with customization options and ordered samples.
    """
    gene_data = summary_df[summary_df['Target Name'] == gene].copy()
    
    if gene_data.empty:
        return None

    # Apply custom sample order
    gene_data['Sample Name'] = pd.Categorical(gene_data['Sample Name'], categories=sample_order, ordered=True)
    gene_data = gene_data.sort_values('Sample Name')

    # Calculate the error bar extent (upper and lower bounds)
    gene_data['Error_Lower'] = gene_data['Relative Expression (2^-ddCt)'] - gene_data['Fold Change Min']
    gene_data['Error_Upper'] = gene_data['Fold Change Max'] - gene_data['Relative Expression (2^-ddCt)']
    
    # Combine upper and lower into one array for matplotlib errorbar
    errors = gene_data[['Error_Lower', 'Error_Upper']].values.T if plot_options['show_error_bars'] else None

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create the bar plot with error bars
    ax.bar(
        gene_data['Sample Name'], 
        gene_data['Relative Expression (2^-ddCt)'], 
        yerr=errors, 
        capsize=plot_options['error_capsize'], 
        color=plot_options['bar_color'], 
        alpha=0.7
    )
    
    # Draw a line at y=1 (Reference/Control Expression Level)
    ax.axhline(1, color='red', linestyle='--', linewidth=1, label='Reference Level (1.0)')

    ax.set_title(f'{plot_options["title_prefix"]} {gene}', fontsize=16, fontweight='bold')
    ax.set_ylabel('Relative Expression (Fold Change)', fontsize=12)
    ax.set_xlabel('Sample Group', fontsize=12)
    
    # Apply custom Y-axis limit
    if plot_options['y_max'] > 1.0:
        ax.set_ylim(0, plot_options['y_max'])
    
    # Style improvements
    ax.grid(axis='y', linestyle=':', alpha=0.6)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    return fig

# ------------------------------------------------------------
# MAIN STREAMLIT APP
# ------------------------------------------------------------

def main():
    st.title("🧬 Interactive qPCR Analysis Tool")
    st.markdown("Upload multiple raw qPCR results files (Excel or CSV) for combined 2^-ΔΔCt analysis.")
    
    # --- 1. Multi-File Upload ---
    uploaded_files = st.file_uploader("Choose one or more files (Excel or CSV)", 
                                      type=['xlsx', 'xls', 'csv'],
                                      accept_multiple_files=True)
    
    if 'data_state' not in st.session_state:
        st.session_state.data_state = None

    if uploaded_files:
        
        # --- 2. Data Loading and Consolidation ---
        all_combined_data_frames = []
        for uploaded_file in uploaded_files:
            file_extension = uploaded_file.name.split('.')[-1].lower()
            with st.spinner(f"Parsing file: {uploaded_file.name}..."):
                standardized_df = parse_qpcr_data(uploaded_file, file_extension) 
            
            if standardized_df is not None and not standardized_df.empty:
                all_combined_data_frames.append(standardized_df)

        if not all_combined_data_frames:
            st.error("No valid data could be extracted from any uploaded file.")
            st.session_state.data_state = None
            return

        all_data = pd.concat(all_combined_data_frames, ignore_index=True)
        standardized_data = standardize_sample_names(all_data)
        cleaned_data = clean_data(standardized_data)
        
        if cleaned_data.empty:
            st.error("No valid Ct data remaining after cleaning.")
            st.session_state.data_state = None
            return
            
        st.success(f"Successfully combined and cleaned data from {len(all_combined_data_frames)} file(s).")
        
        st.session_state.data_state = cleaned_data
        
        # ------------------------------------------------------------------------------------------------
        # --- 3. INTERACTIVE DATA EDITING (NEW SECTION) ---
        # ------------------------------------------------------------------------------------------------
        
        st.subheader("Interactive Sample Editor")
        st.markdown("Edit the name and order of your sample groups before analysis.")

        original_samples = sorted(cleaned_data['Sample Name'].unique().tolist())
        
        # Initialize map in session state
        if 'sample_map' not in st.session_state or st.session_state.original_samples != original_samples:
             st.session_state.sample_map = pd.DataFrame({
                'Original Name': original_samples,
                'New Name': original_samples
            })
             # Store current samples to check if we need to reset the map
             st.session_state.original_samples = original_samples

        col1, col2 = st.columns([1, 1])
        with col1:
            st.markdown("##### 1. Rename Sample Groups")
            
            # Data Editor for Renaming
            edited_map = st.data_editor(
                st.session_state.sample_map,
                column_config={
                    "Original Name": st.column_config.TextColumn("Original Name", disabled=True),
                    "New Name": st.column_config.TextColumn("New Name", required=True)
                },
                use_container_width=True,
                num_rows="fixed",
                key="sample_renamer"
            )
            
            # Update the sample map in session state
            st.session_state.sample_map = edited_map
            
            # Apply renaming to the main data
            renaming_dict = pd.Series(edited_map['New Name'].values, index=edited_map['Original Name']).to_dict()
            edited_data = cleaned_data.copy()
            edited_data['Sample Name'] = edited_data['Sample Name'].map(renaming_dict)
            
            # Get the new list of unique samples for ordering
            new_samples = sorted(edited_data['Sample Name'].unique().tolist())

        with col2:
            st.markdown("##### 2. Order Sample Groups for Plotting")
            
            # Multi-select for ordering (simplest way to achieve custom order)
            ordered_samples = st.multiselect(
                "Select sample groups in the desired order:",
                options=new_samples,
                default=new_samples,
                help="The order you select here will be the order on the X-axis of the graphs."
            )
            
        if not ordered_samples:
            st.warning("Please select and order your sample groups.")
            return

        # ------------------------------------------------------------------------------------------------
        # --- 4. Parameter Selection and Calculation ---
        # ------------------------------------------------------------------------------------------------

        available_targets = sorted(edited_data['Target Name'].unique().tolist())
        
        st.sidebar.header("Analysis Parameters")

        hk_gene = st.sidebar.selectbox(
            "1. Select Housekeeping Gene (Reference Gene)",
            options=[g for g in available_targets if g != 'NTC'],
            index=0 if len(available_targets) > 0 and available_targets[0] != 'NTC' else 0
        )

        ref_group = st.sidebar.selectbox(
            "2. Select Reference Sample Group (Control)",
            options=ordered_samples,
            index=0
        )
        
        target_genes = [g for g in available_targets if g != hk_gene and g != 'NTC']
        
        st.sidebar.markdown("---")
        st.sidebar.markdown("### 3. Genes to Analyze")
        selected_genes = st.sidebar.multiselect(
            "Select Target Genes (Excluding HK)",
            options=target_genes,
            default=target_genes
        )

        # Run Analysis Button
        if st.sidebar.button("Run ΔΔCt Analysis and Visualize"):
            
            if ref_group not in ordered_samples:
                 st.error("The selected reference group must be included in the ordered samples.")
                 return

            # Filter data to only include the HK gene and the genes selected for analysis
            targets_to_include = selected_genes + [hk_gene]
            analysis_data = edited_data[edited_data['Target Name'].isin(targets_to_include)]
            
            if analysis_data.empty:
                st.warning("Selected genes or housekeeping gene not found in the cleaned data.")
                return
            
            with st.spinner("Performing 2^-ΔΔCt calculation..."):
                summary = calculate_ddct(analysis_data, hk_gene, ref_group)
                st.session_state.summary_data = summary

        if 'summary_data' in st.session_state:
            summary = st.session_state.summary_data
            
            st.subheader("Results: Relative Gene Expression (2^-ΔΔCt)")
            
            # --- Results Display ---
            summary_display = summary.copy()
            summary_display = summary_display.round({
                'dCt_Sample_Mean': 3,
                'dCt_SE': 3,
                'ddCt': 3,
                'Relative Expression (2^-ddCt)': 3,
                'Fold Change Min': 3,
                'Fold Change Max': 3
            })
            
            # Style the reference group
            def highlight_ref(s):
                if s['Sample Name'] == ref_group:
                    return ['background-color: #e6ffe6'] * len(s)
                return [''] * len(s)

            st.dataframe(
                summary_display.style.apply(highlight_ref, axis=1),
                use_container_width=True
            )
            
            # --- Download Final Results (Excel) ---
            st.markdown("##### Download Calculated Results")
            excel_buffer = io.BytesIO()
            with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                summary.to_excel(writer, sheet_name='Relative Expression Summary', index=False)
                # Also save the raw data that was used for the analysis
                edited_data.to_excel(writer, sheet_name='Cleaned Raw Data', index=False)
            
            excel_buffer.seek(0)
            
            st.download_button(
                label="⬇️ Download Results & Graph Data (Excel)",
                data=excel_buffer.getvalue(),
                file_name=f"qpcr_final_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            st.markdown("---")
            
            # ------------------------------------------------------------------------------------------------
            # --- 5. VISUALIZATION & CUSTOMIZATION (NEW SECTION) ---
            # ------------------------------------------------------------------------------------------------
            
            st.subheader("Customizable Relative Expression Graphs")
            
            if selected_genes:
                
                # Setup session state for per-gene customization
                if 'plot_options' not in st.session_state:
                    st.session_state.plot_options = {}

                with st.container():
                    st.markdown("Adjust plotting options below, then click **Run Analysis** again to re-render.")
                    
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    zip_buffer = io.BytesIO()
                    
                    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                        
                        for gene in selected_genes:
                            
                            # Initialize default options for a new gene
                            if gene not in st.session_state.plot_options:
                                st.session_state.plot_options[gene] = {
                                    'bar_color': '#1f77b4',
                                    'y_max': round(summary[summary['Target Name'] == gene]['Fold Change Max'].max() * 1.2, 1) or 5.0,
                                    'show_error_bars': True,
                                    'error_capsize': 5,
                                    'title_prefix': 'Relative Expression of'
                                }

                            with st.expander(f"Customize Plot for Gene: {gene}", expanded=False):
                                
                                # Use columns for cleaner layout
                                col_c1, col_c2, col_c3 = st.columns(3)
                                
                                with col_c1:
                                    st.session_state.plot_options[gene]['bar_color'] = st.color_picker(
                                        'Bar Color', 
                                        st.session_state.plot_options[gene]['bar_color'],
                                        key=f"color_{gene}"
                                    )
                                with col_c2:
                                    max_fold_change = summary[summary['Target Name'] == gene]['Fold Change Max'].max()
                                    default_max = round(max_fold_change * 1.2, 1) if pd.notna(max_fold_change) else 5.0
                                    
                                    st.session_state.plot_options[gene]['y_max'] = st.slider(
                                        'Y-Axis Max Limit (0-based)',
                                        min_value=1.0, 
                                        max_value=max(default_max * 2, 10.0), 
                                        value=st.session_state.plot_options[gene]['y_max'], 
                                        step=0.1,
                                        key=f"ymax_{gene}"
                                    )
                                with col_c3:
                                    st.session_state.plot_options[gene]['show_error_bars'] = st.checkbox(
                                        'Show Error Bars (Min/Max FC)',
                                        value=st.session_state.plot_options[gene]['show_error_bars'],
                                        key=f"error_{gene}"
                                    )
                                    st.session_state.plot_options[gene]['title_prefix'] = st.text_input(
                                        'Plot Title Prefix',
                                        st.session_state.plot_options[gene]['title_prefix'],
                                        key=f"prefix_{gene}"
                                    )

                            # Generate the graph with current options
                            fig = create_custom_bar_graph(summary, gene, ordered_samples, st.session_state.plot_options[gene])
                            
                            if fig:
                                st.pyplot(fig)
                                
                                # Save figure to zip for batch download
                                img_buffer = io.BytesIO()
                                fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                                img_buffer.seek(0)
                                zip_file.writestr(f"{gene}_expression.png", img_buffer.getvalue())
                                plt.close(fig) # Close figure to free memory
                                
                    zip_buffer.seek(0)
                    
                    st.download_button(
                        label="📥 Download All Graphs (ZIP)",
                        data=zip_buffer.getvalue(),
                        file_name=f"qpcr_graphs_{timestamp}.zip",
                        mime="application/zip"
                    )
            else:
                st.warning("Please select at least one gene in the sidebar to generate graphs.")

if __name__ == "__main__":
    main()
