"""
qpcr_analysis_app.py — Dynamic and Interactive qPCR Analysis Web Application
-----------------------------------------------------------------------------
Features:
  - User-driven file upload for Excel (.xlsx) or CSV (.csv) files.
  - Interactive configuration of Sample Group and Target Gene order using native Streamlit controls (buttons/select).
  - Flexible parsing: scans data frames for Sample Name, Target Name, Ct columns.
  - Calculates relative gene expression (2^-ΔΔCt) with robust statistics.
  - Generates downloadable, publication-ready results and bar graphs with custom ordering.

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
    page_title="Interactive qPCR Analysis Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ------------------------------------------------------------
# NATIVE STREAMLIT ORDERING UTILITIES (REPLACEMENT FOR ST_ANT)
# ------------------------------------------------------------

def manage_order_state(key, items_list):
    """Initializes or resets the order in session state based on current data."""
    if key not in st.session_state or set(st.session_state[key]) != set(items_list):
        st.session_state[key] = items_list

def interactive_reorder_ui(title, key, all_items):
    """
    Renders a native Streamlit UI for reordering items.
    Updates st.session_state[key] directly.
    """
    
    st.sidebar.markdown(f"**{title}**")
    
    current_order = st.session_state.get(key, all_items)
    
    if not current_order:
        st.sidebar.caption("No items to order.")
        return

    # 1. Select item to move
    col1, col2 = st.sidebar.columns([0.7, 0.3])
    
    item_to_move = col1.selectbox(
        "Select item to move:",
        options=current_order,
        label_visibility="collapsed",
        key=f'{key}_select'
    )

    # Find the index of the selected item
    try:
        idx = current_order.index(item_to_move)
    except ValueError:
        # Should not happen if item is in current_order, but handles safety
        return

    # 2. Movement buttons
    with col2:
        st.write("") # Spacer for alignment
        up_btn = st.button("⬆️", key=f'{key}_up', use_container_width=True)
        down_btn = st.button("⬇️", key=f'{key}_down', use_container_width=True)

    
    # 3. Handle movement logic
    new_order = list(current_order)
    moved = False

    if up_btn and idx > 0:
        # Swap with the item above
        new_order[idx], new_order[idx - 1] = new_order[idx - 1], new_order[idx]
        moved = True
    elif down_btn and idx < len(new_order) - 1:
        # Swap with the item below
        new_order[idx], new_order[idx + 1] = new_order[idx + 1], new_order[idx]
        moved = True

    if moved:
        # Update the session state immediately
        st.session_state[key] = new_order
        # Need to rerun to update the UI with the new order in the selectbox
        st.rerun()

    # 4. Display current order
    st.sidebar.markdown(f"<p style='font-size:12px; margin-top:5px; margin-bottom:0;'>Current Order:</p>", unsafe_allow_html=True)
    st.sidebar.markdown(f"<p style='font-size:14px; padding: 5px; border: 1px solid #ddd; border-radius: 4px;'>{ ' → '.join(current_order) }</p>", unsafe_allow_html=True)


# ------------------------------------------------------------
# FLEXIBLE DATA SCANNING FUNCTIONS
# ------------------------------------------------------------

@st.cache_data
def scan_for_columns(df):
    """Scan DataFrame for Sample Name, Target Name, and Ct columns."""
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
    """Loads data from uploaded file objects, combines, cleans, and prepares."""
    st.info("Starting data parsing and combination...")
    all_data = []

    for i, uploaded_file in enumerate(uploaded_files):
        file_name = uploaded_file.name
        header_row_index = header_rows[i] - 1
        
        try:
            if file_name.endswith('.xlsx'):
                df = pd.read_excel(uploaded_file, header=header_row_index, engine='openpyxl')
            elif file_name.endswith('.csv'):
                df = pd.read_csv(uploaded_file, header=header_row_index)
            else:
                st.warning(f"Skipping file **{file_name}**: Unsupported format.")
                continue

            df.columns = [col.strip() for col in df.columns]
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
                    st.success(f"✅ Processed file: **{file_name}**")
                else:
                    st.warning(f"⚠️ File **{file_name}** had no valid Ct data after cleaning. Skipped.")
            else:
                st.error(f"❌ Failed to find mandatory columns (Sample Name, Target Name, Ct) in file **{file_name}**.")

        except Exception as e:
            st.error(f"Error loading file **{file_name}**: {e}")

    if not all_data:
        st.error("No valid qPCR data could be loaded for analysis.")
        return None, None, None, None

    combined_df = pd.concat(all_data, ignore_index=True)
    combined_df['Group'] = combined_df['Sample Name'].apply(lambda x: str(x).strip())
    
    samples = combined_df['Sample Name'].unique().tolist()
    targets = combined_df['Target Name'].unique().tolist()
    groups = combined_df['Group'].unique().tolist()
    
    st.success(f"Data loading complete. Found {len(samples)} unique samples and {len(targets)} unique targets.")
    
    return combined_df, samples, targets, groups

# ------------------------------------------------------------
# CALCULATION & SORTING FUNCTIONS
# ------------------------------------------------------------

def calculate_ddct(data, control_gene, reference_group):
    """Performs the 2^-ΔΔCt calculation and T-Test for significance."""
    
    # Check for selected gene/group presence before calculation
    if control_gene not in data['Target Name'].unique():
        st.error(f"Control gene **'{control_gene}'** not found in data.")
        return None
    if reference_group not in data['Group'].unique():
        st.error(f"Reference group **'{reference_group}'** not found in data.")
        return None

    # 1. Group by Sample and Target, calculate Mean Ct
    ct_mean = data.groupby(['Group', 'Target Name'])['Ct'].mean().reset_index().rename(columns={'Ct': 'Mean Ct'})
    
    # 2. Calculate ΔCt
    control_cts = ct_mean[ct_mean['Target Name'] == control_gene].set_index('Group')['Mean Ct']
    if control_cts.empty: return None

    ct_mean['Control Ct'] = ct_mean['Group'].map(control_cts)
    ct_mean = ct_mean.dropna(subset=['Control Ct']).copy() 
    ct_mean['ΔCt'] = ct_mean['Mean Ct'] - ct_mean['Control Ct']

    # 3. Calculate ΔΔCt
    ref_dct = ct_mean[
        (ct_mean['Group'] == reference_group) & (ct_mean['Target Name'] != control_gene)
    ].set_index('Target Name')['ΔCt']
    if ref_dct.empty:
        st.error(f"The reference group **'{reference_group}'** did not have detectable ΔCt values for any target gene.")
        return None
        
    ct_mean['Reference ΔCt'] = ct_mean.apply(lambda row: ref_dct.get(row['Target Name']), axis=1)
    summary_df = ct_mean[ct_mean['Target Name'] != control_gene].copy()
    summary_df = summary_df.dropna(subset=['Reference ΔCt']).copy()
    summary_df['ΔΔCt'] = summary_df['ΔCt'] - summary_df['Reference ΔCt']
    
    # 4. Fold Change
    summary_df['Fold Change (2^-ΔΔCt)'] = 2**(-summary_df['ΔΔCt'])
    
    # 5. Calculate Standard Error (SE) for ΔΔCt
    raw_data_agg = data.groupby(['Group', 'Target Name'])['Ct'].agg(['mean', 'std', 'count']).reset_index().rename(columns={'mean': 'Mean Ct', 'std': 'Std Ct', 'count': 'N'})
    control_agg = raw_data_agg[raw_data_agg['Target Name'] == control_gene].set_index('Group')[['Std Ct', 'N']].rename(columns={'Std Ct': 'Std Ct Control', 'N': 'N Control'})
    full_agg = raw_data_agg.merge(control_agg, on='Group', how='left')
    full_agg = full_agg[full_agg['Target Name'] != control_gene].dropna(subset=['Std Ct Control']).copy()

    # SE(ΔCt)
    full_agg['SE ΔCt'] = np.sqrt((full_agg['Std Ct']**2 / full_agg['N']) + (full_agg['Std Ct Control']**2 / full_agg['N Control']))
    
    ref_se_dct = full_agg[full_agg['Group'] == reference_group].set_index('Target Name')['SE ΔCt']
    full_agg['Ref SE ΔCt'] = full_agg['Target Name'].map(ref_se_dct)
    full_agg = full_agg.dropna(subset=['Ref SE ΔCt']).copy()
    
    # SE ΔΔCt
    full_agg['SE ΔΔCt'] = np.sqrt(full_agg['SE ΔCt']**2 + full_agg['Ref SE ΔCt']**2)
    
    summary_df_se = summary_df.merge(full_agg[['Group', 'Target Name', 'SE ΔΔCt']], on=['Group', 'Target Name'], how='left')
    summary_df_se = summary_df_se.set_index(['Group', 'Target Name'])
    
    summary_df_se['Upper Bound'] = 2**(-(summary_df_se['ΔΔCt'] - summary_df_se['SE ΔΔCt']))
    summary_df_se['Lower Bound'] = 2**(-(summary_df_se['ΔΔCt'] + summary_df_se['SE ΔΔCt']))
    
    # 6. T-Test
    p_values = {}
    for (group, target), _ in summary_df_se.iterrows():
        if group == reference_group:
            p_values[(group, target)] = np.nan
            continue
            
        # Helper to get ΔCt for each biological replicate (Sample Name)
        def get_raw_dct_for_group(data_df, group_name, target_name, control_name):
            raw_dcts = []
            for sample in data_df[data_df['Group'] == group_name]['Sample Name'].unique():
                target_ct = data_df[(data_df['Sample Name'] == sample) & (data_df['Target Name'] == target_name)]['Ct'].mean()
                control_ct = data_df[(data_df['Sample Name'] == sample) & (data_df['Target Name'] == control_name)]['Ct'].mean()
                if not pd.isna(target_ct) and not pd.isna(control_ct):
                    raw_dcts.append(target_ct - control_ct)
            return raw_dcts

        group_dct_raw = get_raw_dct_for_group(data, group, target, control_gene)
        ref_dct_raw = get_raw_dct_for_group(data, reference_group, target, control_gene)

        if len(group_dct_raw) > 1 and len(ref_dct_raw) > 1:
            _, p_val = stats.ttest_ind(group_dct_raw, ref_dct_raw, equal_var=False)
            p_values[(group, target)] = p_val
        else:
            p_values[(group, target)] = np.nan
            
    summary_df_se['P-Value'] = pd.Series(p_values)
    
    # 7. Final Formatting
    final_summary = summary_df_se.reset_index()
    final_summary['Significant'] = final_summary['P-Value'].apply(
        lambda x: 'Yes (p < 0.05)' if x < 0.05 else ('No' if not pd.isna(x) else 'N/A')
    )
    
    return final_summary.filter(items=[
        'Group', 'Target Name', 'Mean Ct', 'ΔCt', 'ΔΔCt', 
        'Fold Change (2^-ΔΔCt)', 'SE ΔΔCt', 'Upper Bound', 'Lower Bound', 'P-Value', 'Significant'
    ], axis=1)

def apply_user_ordering(df, target_order, group_order):
    """Applies user-defined target and group order, with numerical fallback for groups."""
    
    if df is None or df.empty:
        return pd.DataFrame()

    # 1. Apply Target Gene Order (Categorical)
    known_targets = df['Target Name'].unique().tolist()
    full_target_order = [t for t in target_order if t in known_targets] + [t for t in known_targets if t not in target_order]
    
    df['Target_Order'] = pd.Categorical(df['Target Name'], categories=full_target_order, ordered=True)
    
    # 2. Apply Group Order (Categorical for Plotting)
    known_groups = df['Group'].unique().tolist()
    full_group_order = [g for g in group_order if g in known_groups] + [g for g in known_groups if g not in group_order]
    
    df['Group_Order'] = pd.Categorical(df['Group'], categories=full_group_order, ordered=True)
    
    # 3. Create a Numerical Sort Key for Groups
    def numerical_sort_key(group_name):
        # Extracts numerical parts and converts them to integers for correct numerical sorting
        parts = re.split('(\d+)', str(group_name))
        return [int(text) if text.isdigit() else text.lower() for text in parts]
    
    df['Numerical_Key'] = df['Group'].apply(numerical_sort_key)
    
    # 4. Final Sorting
    # Primary sort: User-defined Target Order (Groups the genes)
    # Secondary sort: Numerical_Key (Ensures Sample 1 comes before Sample 10)
    # Tertiary sort: User-defined Group_Order (Ensures the categorical plot order is maintained for stability)
    df_sorted = df.sort_values(
        by=['Target_Order', 'Numerical_Key'],
        kind='mergesort' 
    ).drop(columns=['Target_Order', 'Group_Order', 'Numerical_Key']).reset_index(drop=True)
    
    return df_sorted

# ------------------------------------------------------------
# PLOTTING FUNCTIONS 
# ------------------------------------------------------------

def create_bar_graph(summary, gene, group_order):
    """
    Creates a bar plot for the specified gene, respecting the user's group order.
    """
    plot_data = summary[summary['Target Name'] == gene].copy()
    
    if plot_data.empty:
        return None
    
    # Apply group order to the plotting data
    plot_data['Group'] = pd.Categorical(plot_data['Group'], categories=group_order, ordered=True)
    plot_data = plot_data.sort_values('Group')
        
    ref_group = plot_data[plot_data['Fold Change (2^-ΔΔCt)'] == 1]['Group'].iloc[0] if not plot_data[plot_data['Fold Change (2^-ΔΔCt)'] == 1].empty else None
    
    plot_data['Error'] = plot_data['Upper Bound'] - plot_data['Fold Change (2^-ΔΔCt)']
    if ref_group in plot_data['Group'].values:
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
            # Adjust position slightly above the upper bound
            ax.text(i, row['Upper Bound'] + 0.1, '*', ha='center', va='bottom', fontsize=18, color='darkred', weight='bold')
        
    ax.axhline(1, color='gray', linestyle='--', linewidth=1, label='Reference (Fold Change = 1)')
    ax.set_ylabel('Fold Change (2⁻ΔΔCt) Relative to Reference')
    ax.set_xlabel('Sample Group')
    ax.set_title(f'Relative Gene Expression: {gene}', fontsize=16)
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
    
    st.title("🧬 Interactive qPCR 2⁻ΔΔCt Analysis Tool")
    
    # --- SIDEBAR CONFIGURATION ---
    
    st.sidebar.header("1. Upload Data")
    uploaded_files = st.sidebar.file_uploader(
        "Upload one or more qPCR results files",
        type=['csv', 'xlsx'],
        accept_multiple_files=True
    )
    
    if not uploaded_files:
        st.info("Please upload your raw data files in the sidebar to proceed.")
        return

    # 1.1 Header Row Selector
    st.sidebar.markdown("---")
    st.sidebar.subheader("File Configuration")
    header_rows = []
    
    for i, file in enumerate(uploaded_files):
        header_row = st.sidebar.number_input(
            f"Header Row (1-based) for: **{file.name}**",
            min_value=1, 
            value=1,
            step=1, 
            key=f'header_{i}'
        )
        header_rows.append(header_row)

    # 2. Load Data (Caches the result)
    data, samples, targets, groups = load_and_preprocess_data(uploaded_files, header_rows)
    
    if data is None:
        return 
        
    st.sidebar.markdown("---")
    st.sidebar.header("2. Set Analysis Parameters")
    
    # 2.1 Select Housekeeping/Control Gene
    default_control_index = targets.index('ACTIN') if 'ACTIN' in targets else (0 if targets else None)
    control_gene = st.sidebar.selectbox(
        "Select Housekeeping/Control Gene:",
        targets,
        index=default_control_index if default_control_index is not None else 0
    )
    
    # 2.2 Select Reference Group
    reference_group = st.sidebar.selectbox(
        "Select Reference Group (Fold Change = 1):",
        groups,
        index=0
    )
    
    st.sidebar.markdown("---")
    st.sidebar.header("3. Interactive Ordering")
    
    # --- Interactive Gene Ordering ---
    target_genes = [t for t in targets if t != control_gene]
    
    # Ensure state is initialized
    manage_order_state('target_order', target_genes)
    manage_order_state('group_order', groups)
    
    # Use native Streamlit UI for reordering
    with st.sidebar.expander("Configure Target Gene Order", expanded=True):
        interactive_reorder_ui("Target Gene Order", 'target_order', target_genes)
    
    with st.sidebar.expander("Configure Sample Group Order", expanded=True):
        interactive_reorder_ui("Sample Group Order", 'group_order', groups)
    
    st.sidebar.markdown("---")
    
    # 2.3 Calculate Button
    if st.sidebar.button("🔬 Start 2⁻ΔΔCt Calculation", use_container_width=True):
        
        # Perform calculation
        summary = calculate_ddct(data, control_gene, reference_group)
        
        if summary is not None and not summary.empty:
            
            # Apply user sorting to the summary table
            final_summary_sorted = apply_user_ordering(
                summary, 
                st.session_state['target_order'], 
                st.session_state['group_order']
            )
            
            st.session_state['final_summary'] = final_summary_sorted
            st.session_state['targets'] = targets
            st.session_state['control_gene'] = control_gene
            st.success("Calculation and ordering complete! Results are below.")
        
    # --- Display Results ---
    if 'final_summary' in st.session_state:
        final_summary_sorted = st.session_state['final_summary']
        control_gene = st.session_state['control_gene']

        st.subheader("4. Results: Custom Sorted Relative Gene Expression Summary")
        
        # Display summary table with formatting
        display_summary = final_summary_sorted.copy()
        numeric_cols = ['Mean Ct', 'ΔCt', 'ΔΔCt', 'Fold Change (2^-ΔΔCt)', 'SE ΔΔCt', 'Upper Bound', 'Lower Bound', 'P-Value']
        for col in numeric_cols:
            if col in display_summary.columns:
                display_summary[col] = display_summary[col].round(3)
        
        st.markdown(f"**Target Genes** are grouped and ordered by your selection. Within each gene, **Groups** are sorted numerically for clean viewing.")
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
            file_name=f"qpcr_summary_sorted_{timestamp}.csv",
            mime="text/csv"
        )
        
        # --- Graphing Section ---
        st.markdown("---")
        st.subheader("5. Visualization & Graphs")
        
        # Use the ordered targets for the multiselect
        target_genes_ordered = [t for t in st.session_state['target_order'] if t != control_gene]
        
        selected_genes = st.multiselect(
            "Select target genes to plot (Will be plotted in the selected order):",
            options=target_genes_ordered,
            default=target_genes_ordered,
            key='plot_genes_select'
        )
        
        if st.button("📈 Generate Graphs", use_container_width=True):
            if selected_genes:
                with st.spinner("Creating graphs..."):
                    
                    zip_buffer = io.BytesIO()
                    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                        
                        current_group_order = st.session_state['group_order']

                        for gene in selected_genes:
                            fig = create_bar_graph(final_summary_sorted, gene, current_group_order)
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


if __name__ == "__main__":
    main()
