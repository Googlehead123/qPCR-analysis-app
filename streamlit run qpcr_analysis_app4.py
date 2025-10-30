"""
qpcr_analysis_app.py — Interactive and Customizable qPCR Analysis Tool
---------------------------------------------------------------------
Features:
  - Multi-file parsing: scans all sheets and multiple uploaded Excel/CSV files.
  - Robust column matching (Sample Name, Target Name, Ct).
  - Standardizes numeric sample names (e.g., '1' and '1.0' become '1').
  - INTERACTIVE DATA EDITING: Allows renaming, reordering, and filtering of sample groups.
  - Calculates relative gene expression (2^-ΔΔCt) grouped by target gene.
  - Customizable bar graphs (color, Y-axis limits, titles) for each gene using Plotly.
  - Downloadable final results and cleaned raw data as an Excel file.

Usage:
    streamlit run qpcr_analysis_app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import io
from datetime import datetime
import zipfile
import re
import plotly.express as px
import plotly.graph_objects as go

# Set page config
st.set_page_config(
    page_title="Customizable qPCR Analysis Tool",
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
    
    # Reset file pointer
    uploaded_file.seek(0) 

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
            # Need to find header in CSV similarly to Excel
            df = pd.read_csv(uploaded_file, header=None)
            uploaded_file.seek(0) # Reset pointer for second read
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
            # Select and copy only the required columns
            standardized_data = standardized_data[['Sample Name', 'Target Name', 'Ct']].copy()
            return standardized_data
        
        return None

    except Exception as e:
        st.error(f"An error occurred during parsing {uploaded_file.name}: {e}")
        return None

# ------------------------------------------------------------
# DATA CLEANING AND STANDARDIZATION
# ------------------------------------------------------------

def standardize_sample_names(df):
    """
    Standardizes numeric sample names (e.g., convert 1.0 to 1).
    """
    sample_col = 'Sample Name'
    # Attempt to convert to string and strip whitespace
    s = df[sample_col].astype(str).str.strip().copy()
    numeric_s = pd.to_numeric(s, errors='coerce')
    is_numeric = numeric_s.notna()
    
    try:
        # Convert valid numbers to integer strings (e.g., '1.0' -> '1')
        df.loc[is_numeric, sample_col] = numeric_s[is_numeric].round(0).astype(int).astype(str)
    except Exception:
        # Fallback for complex types, just convert to string and remove trailing .0
        df.loc[is_numeric, sample_col] = numeric_s[is_numeric].astype(str).str.replace(r'\.0$', '', regex=True)
    
    df[sample_col] = df[sample_col].astype(str).str.strip()
    return df

@st.cache_data
def clean_data(df):
    """
    Filters out non-numeric Ct values, ensures columns are correct types, 
    and filters out high Ct values.
    """
    df['Ct'] = pd.to_numeric(df['Ct'], errors='coerce')
    # Filter out NaN Ct values and values above 35 (considered indeterminate/non-detected)
    df = df.dropna(subset=['Ct'])
    df = df[df['Ct'] <= 35] 
    
    # Ensure name columns are clean strings
    df['Target Name'] = df['Target Name'].astype(str).str.strip()
    df['Sample Name'] = df['Sample Name'].astype(str).str.strip()
    
    # Filter out empty names
    df = df[df['Sample Name'].str.len() > 0]
    df = df[df['Target Name'].str.len() > 0]
    
    return df

# ------------------------------------------------------------
# ΔΔCt CALCULATION FUNCTION (CACHE FOR PERFORMANCE)
# ------------------------------------------------------------

@st.cache_data
def calculate_ddct(df, hk_gene, ref_group):
    """
    Performs the 2^-ΔΔCt calculation using the chosen housekeeping gene and reference group.
    """
    
    # 1. Calculate ΔCt (Ct_Target - Ct_HK)
    # Calculate Mean Ct for the Housekeeping (HK) gene per Sample
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
    
    # Filter out the HK gene from the dCt summary for the next steps
    dCt_summary = dCt_summary[dCt_summary['Target Name'] != hk_gene].copy()
    
    # 2. Calculate ΔΔCt (ΔCt_Sample - ΔCt_Reference)
    
    # Get the ΔCt_Ref_Mean for each target gene from the reference group
    ref_dCt = dCt_summary[dCt_summary['Sample Name'] == ref_group][['Target Name', 'dCt_Mean']]
    ref_dCt = ref_dCt.rename(columns={'dCt_Mean': 'dCt_Ref_Mean'})
    
    summary = pd.merge(dCt_summary, ref_dCt, on='Target Name', how='left')
    
    # Calculate ΔΔCt (Handle potential NaNs if a gene is missing in the reference group)
    summary['ddCt'] = summary['dCt_Sample_Mean'] - summary['dCt_Ref_Mean']
    
    # 3. Calculate Relative Expression (2^-ΔΔCt)
    summary['Relative Expression (2^-ddCt)'] = 2**(-summary['ddCt'])
    
    # Calculate Fold Change Min and Max based on SE
    summary['Fold Change Min'] = 2**(-(summary['ddCt'] + summary['dCt_SE']))
    summary['Fold Change Max'] = 2**(-(summary['ddCt'] - summary['dCt_SE']))
    
    final_summary = summary.drop(columns=['dCt_Ref_Mean']).rename(columns={'dCt_Sample_Mean': 'dCt_Mean'})
    
    # For the reference group, ensure fold change is exactly 1.0 where possible
    final_summary.loc[final_summary['Sample Name'] == ref_group, 'Relative Expression (2^-ddCt)'] = 1.0
    final_summary.loc[final_summary['Sample Name'] == ref_group, 'ddCt'] = 0.0
    
    return final_summary

# ------------------------------------------------------------
# INTERACTIVE GRAPHING FUNCTION (USING PLOTLY)
# ------------------------------------------------------------

def create_custom_bar_graph(summary_df, gene, sample_order, plot_options, ref_group):
    """
    Creates a Plotly bar graph for a single gene with customization options and ordered samples.
    """
    gene_data = summary_df[summary_df['Target Name'] == gene].copy()
    
    if gene_data.empty:
        return None

    # Apply custom sample order
    gene_data['Sample Name'] = pd.Categorical(gene_data['Sample Name'], categories=sample_order, ordered=True)
    gene_data = gene_data.sort_values('Sample Name')

    # Prepare error bars
    error_y = None
    if plot_options['show_error_bars']:
        gene_data['Error_Lower'] = gene_data['Relative Expression (2^-ddCt)'] - gene_data['Fold Change Min']
        gene_data['Error_Upper'] = gene_data['Fold Change Max'] - gene_data['Relative Expression (2^-ddCt)']
        
        error_y = {
            'type': 'data',
            'symmetric': False,
            'array': gene_data['Error_Upper'].tolist(),
            'arrayminus': gene_data['Error_Lower'].tolist()
        }

    # Create the Plotly figure
    fig = go.Figure(data=[
        go.Bar(
            x=gene_data['Sample Name'],
            y=gene_data['Relative Expression (2^-ddCt)'],
            error_y=error_y,
            marker_color=plot_options['bar_color'],
            name='Fold Change',
            # Add tooltip data
            hovertemplate=(
                '<b>Sample</b>: %{x}<br>' +
                '<b>Fold Change</b>: %{y:.3f}<br>' +
                '<b>dCt SE</b>: %{customdata[0]:.3f}<br>' +
                '<extra></extra>' # Hides the trace name
            ),
            customdata=gene_data[['dCt_SE']].values
        )
    ])
    
    # Add Reference Line at y=1
    fig.add_hline(
        y=1, 
        line_dash="dot", 
        line_color="red", 
        annotation_text=f"Reference ({ref_group})", 
        annotation_position="bottom right"
    )

    # Apply layout customizations
    fig.update_layout(
        title={
            'text': f"{plot_options['title_text']}", 
            'y':0.9, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'
        },
        yaxis_title=plot_options['y_label'],
        xaxis_title=plot_options['x_label'],
        # Apply custom Y-axis limit
        yaxis_range=[0, plot_options['y_max']],
        height=500
    )
    
    return fig

# ------------------------------------------------------------
# MAIN STREAMLIT APP
# ------------------------------------------------------------

def main():
    st.title("🧬 Interactive Customizable qPCR Analysis Tool")
    st.markdown("Upload multiple raw qPCR results files (Excel or CSV) for combined 2^-ΔΔCt analysis.")
    
    # --- 1. Multi-File Upload ---
    uploaded_files = st.file_uploader("Choose one or more files (Excel or CSV)", 
                                      type=['xlsx', 'xls', 'csv'],
                                      accept_multiple_files=True)
    
    # Initialize session state for data flow
    if 'cleaned_data' not in st.session_state: st.session_state.cleaned_data = pd.DataFrame()
    if 'edited_data' not in st.session_state: st.session_state.edited_data = pd.DataFrame()
    if 'summary_data' not in st.session_state: st.session_state.summary_data = pd.DataFrame()

    if uploaded_files:
        
        # --- Data Loading and Consolidation ---
        all_combined_data_frames = []
        for uploaded_file in uploaded_files:
            file_extension = uploaded_file.name.split('.')[-1].lower()
            with st.spinner(f"Parsing file: {uploaded_file.name}..."):
                standardized_df = parse_qpcr_data(uploaded_file, file_extension) 
            
            if standardized_df is not None and not standardized_df.empty:
                all_combined_data_frames.append(standardized_df)

        if not all_combined_data_frames:
            st.error("No valid data could be extracted from any uploaded file.")
            st.session_state.cleaned_data = pd.DataFrame()
            return

        all_data = pd.concat(all_combined_data_frames, ignore_index=True)
        standardized_data = standardize_sample_names(all_data)
        cleaned_data = clean_data(standardized_data)
        
        if cleaned_data.empty:
            st.error("No valid Ct data remaining after cleaning.")
            st.session_state.cleaned_data = pd.DataFrame()
            return
            
        st.session_state.cleaned_data = cleaned_data
        st.success(f"Successfully combined and cleaned raw data from {len(all_combined_data_frames)} file(s).")
        
        # Download Cleaned Raw Data
        with st.expander("Download Cleaned Raw Data", expanded=False):
            st.markdown("This is the data after initial parsing, standardization, and filtering (Ct > 35 removed).")
            excel_buffer_clean = io.BytesIO()
            with pd.ExcelWriter(excel_buffer_clean, engine='xlsxwriter') as writer:
                cleaned_data.to_excel(writer, sheet_name='Cleaned Raw Data', index=False)
            excel_buffer_clean.seek(0)
            st.download_button(
                label="⬇️ Download Cleaned Raw Data (Excel)",
                data=excel_buffer_clean.getvalue(),
                file_name=f"qpcr_cleaned_raw_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            st.dataframe(cleaned_data.head())


        st.markdown("---")
        # ------------------------------------------------------------------------------------------------
        # --- 2. INTERACTIVE DATA EDITING (Renaming, Filtering, Ordering) ---
        # ------------------------------------------------------------------------------------------------
        
        st.subheader("📝 Interactive Sample Editor")
        st.markdown("Rename, omit, and reorder your sample groups before running the ΔΔCt analysis.")

        original_samples = sorted(cleaned_data['Sample Name'].unique().tolist())
        
        # Initialize map in session state only if data or samples changed
        if 'sample_map' not in st.session_state or st.session_state.get('original_samples_cache') != original_samples:
             st.session_state.sample_map = pd.DataFrame({
                'Original Name': original_samples,
                'New Name': original_samples,
                'Include in Analysis': True # New column for filtering
            })
             st.session_state.original_samples_cache = original_samples # Cache current samples

        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("##### 1. Rename & Omit Sample Groups")
            
            # Data Editor for Renaming and Filtering
            edited_map = st.data_editor(
                st.session_state.sample_map,
                column_config={
                    "Original Name": st.column_config.TextColumn("Original Name", disabled=True),
                    "New Name": st.column_config.TextColumn("New Name", required=True),
                    "Include in Analysis": st.column_config.CheckboxColumn("Include", help="Uncheck to omit sample from analysis")
                },
                use_container_width=True,
                num_rows="fixed",
                key="sample_renamer"
            )
            
            st.session_state.sample_map = edited_map
            
            # Apply renaming and filtering
            filtered_map = edited_map[edited_map['Include in Analysis'] == True]
            renaming_dict = pd.Series(filtered_map['New Name'].values, index=filtered_map['Original Name']).to_dict()
            
            # Filter the main data to only include samples marked for analysis
            samples_to_include = filtered_map['Original Name'].tolist()
            edited_data = cleaned_data[cleaned_data['Sample Name'].isin(samples_to_include)].copy()
            
            # Apply renaming
            edited_data['Sample Name'] = edited_data['Sample Name'].map(renaming_dict)
            st.session_state.edited_data = edited_data

            # Get the new list of unique samples for ordering
            new_samples = sorted(edited_data['Sample Name'].unique().tolist())

        with col2:
            st.markdown("##### 2. Order Sample Groups for Plotting")
            
            # Multi-select for ordering
            ordered_samples = st.multiselect(
                "Select sample groups in the desired order (Drag and drop to reorder):",
                options=new_samples,
                default=new_samples,
                help="The order of samples shown here will be used for all graphs and table presentation."
            )
            st.session_state.ordered_samples = ordered_samples
            
        if not st.session_state.ordered_samples:
            st.warning("Please define the samples to include and their order for analysis.")
            return

        st.markdown("---")
        # ------------------------------------------------------------------------------------------------
        # --- 3. PARAMETER SELECTION AND CALCULATION ---
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
            options=st.session_state.ordered_samples,
            index=st.session_state.ordered_samples.index(st.session_state.ordered_samples[0]) if st.session_state.ordered_samples else 0
        )
        
        target_genes = [g for g in available_targets if g != hk_gene and g != 'NTC']
        
        st.sidebar.markdown("---")
        st.sidebar.markdown("### 3. Genes to Analyze")
        selected_genes = st.sidebar.multiselect(
            "Select Target Genes (Excluding HK)",
            options=target_genes,
            default=target_genes
        )
        st.session_state.selected_genes = selected_genes

        # Run Analysis Button
        if st.sidebar.button("Run ΔΔCt Analysis and Visualize"):
            
            if ref_group not in st.session_state.ordered_samples:
                 st.error("The selected reference group must be included in the ordered samples.")
                 return

            # Filter data to only include the HK gene and the genes selected for analysis
            targets_to_include = selected_genes + [hk_gene]
            analysis_data = st.session_state.edited_data[st.session_state.edited_data['Target Name'].isin(targets_to_include)]
            
            if analysis_data.empty:
                st.warning("Selected genes or housekeeping gene not found in the analyzed data set.")
                return
            
            with st.spinner("Performing 2^-ΔΔCt calculation..."):
                summary = calculate_ddct(analysis_data, hk_gene, ref_group)
                st.session_state.summary_data = summary

        # ------------------------------------------------------------------------------------------------
        # --- 4. RESULTS AND VISUALIZATION ---
        # ------------------------------------------------------------------------------------------------

        if not st.session_state.summary_data.empty:
            summary = st.session_state.summary_data
            
            st.subheader("Final Results: Relative Gene Expression (2^-ΔΔCt)")
            
            # --- Results Display ---
            summary_display = summary.copy()
            
            # Round for display
            summary_display = summary_display.round({
                'dCt_Mean': 3,
                'dCt_SE': 3,
                'ddCt': 3,
                'Relative Expression (2^-ddCt)': 3,
                'Fold Change Min': 3,
                'Fold Change Max': 3
            })
            
            # Sort by Target Gene and then Sample Order
            summary_display['Sample Name'] = pd.Categorical(summary_display['Sample Name'], 
                                                            categories=st.session_state.ordered_samples, 
                                                            ordered=True)
            summary_display = summary_display.sort_values(['Target Name', 'Sample Name'])


            st.dataframe(summary_display, use_container_width=True)
            
            # --- Download Final Results (Excel) ---
            st.markdown("##### Download Final Calculated Results")
            excel_buffer = io.BytesIO()
            with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                summary.to_excel(writer, sheet_name='Relative Expression Summary', index=False)
                st.session_state.edited_data.to_excel(writer, sheet_name='Edited Raw Data', index=False)
            
            excel_buffer.seek(0)
            
            st.download_button(
                label="⬇️ Download Final Analysis & Edited Raw Data (Excel)",
                data=excel_buffer.getvalue(),
                file_name=f"qpcr_final_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            st.markdown("---")
            
            # ------------------------------------------------------------------------------------------------
            # --- 5. VISUALIZATION & CUSTOMIZATION (Per Gene) ---
            # ------------------------------------------------------------------------------------------------
            
            st.subheader("📊 Customizable Gene Expression Graphs")
            
            if st.session_state.selected_genes:
                
                # Setup session state for per-gene customization
                if 'plot_options' not in st.session_state:
                    st.session_state.plot_options = {}

                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                zip_buffer = io.BytesIO()
                
                with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                    
                    for gene in st.session_state.selected_genes:
                        
                        # Initialize default options for a new gene
                        if gene not in st.session_state.plot_options:
                            max_fold_change = summary[summary['Target Name'] == gene]['Fold Change Max'].max()
                            
                            st.session_state.plot_options[gene] = {
                                'bar_color': '#1f77b4',
                                'y_max': round(max_fold_change * 1.2, 1) if pd.notna(max_fold_change) else 5.0,
                                'show_error_bars': True,
                                'title_text': f'Relative Expression of {gene}',
                                'x_label': 'Sample Group',
                                'y_label': 'Relative Expression (Fold Change)'
                            }

                        with st.expander(f"⚙️ Customize Plot for Gene: **{gene}**", expanded=False):
                            
                            # Use columns for cleaner layout
                            col_c1, col_c2, col_c3 = st.columns(3)
                            
                            with col_c1:
                                st.session_state.plot_options[gene]['bar_color'] = st.color_picker(
                                    'Bar Color', 
                                    st.session_state.plot_options[gene]['bar_color'],
                                    key=f"color_{gene}"
                                )
                                st.session_state.plot_options[gene]['show_error_bars'] = st.checkbox(
                                    'Show Error Bars (Min/Max FC)',
                                    value=st.session_state.plot_options[gene]['show_error_bars'],
                                    key=f"error_{gene}"
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
                                st.session_state.plot_options[gene]['x_label'] = st.text_input(
                                    'X-Axis Label',
                                    st.session_state.plot_options[gene]['x_label'],
                                    key=f"xlabel_{gene}"
                                )
                            with col_c3:
                                st.session_state.plot_options[gene]['title_text'] = st.text_input(
                                    'Plot Title',
                                    st.session_state.plot_options[gene]['title_text'],
                                    key=f"title_{gene}"
                                )
                                st.session_state.plot_options[gene]['y_label'] = st.text_input(
                                    'Y-Axis Label',
                                    st.session_state.plot_options[gene]['y_label'],
                                    key=f"ylabel_{gene}"
                                )
                                

                        # Generate the graph with current options
                        fig = create_custom_bar_graph(
                            summary, 
                            gene, 
                            st.session_state.ordered_samples, 
                            st.session_state.plot_options[gene], 
                            ref_group
                        )
                        
                        if fig:
                            st.plotly_chart(fig, use_container_width=True)
                            
                            # Save figure to zip for batch download
                            img_buffer = io.BytesIO()
                            fig.write_image(img_buffer, format='png', scale=2)
                            img_buffer.seek(0)
                            zip_file.writestr(f"{gene}_expression.png", img_buffer.getvalue())
                            
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
