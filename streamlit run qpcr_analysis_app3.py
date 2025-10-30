"""
qpcr_analysis_app.py — Interactive qPCR Analysis Tool
Usage: streamlit run qpcr_analysis_app.py
"""
import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import io
from datetime import datetime
import zipfile

st.set_page_config(page_title="qPCR Analysis", page_icon="🧬", layout="wide")

@st.cache_data
def scan_for_columns(df):
    patterns = {'sample': ['sample name', 'sample', 'samplename', 'name'],
                'target': ['target name', 'target', 'targetname', 'gene', 'assay', 'reporter'],
                'ct': ['ct', 'c t', 'c-t', 'cq', 'value', 'cт', 'c т', 'ст', 'ct mean', 'cт mean']}
    normalize = lambda h: str(h).strip().lower().replace(' ', '').replace('-', '').replace('/', '')
    df_cols = [normalize(col) for col in df.columns]
    found = {}
    for col_type, variations in patterns.items():
        norm_vars = [normalize(v) for v in variations]
        idx = next((i for i, dc in enumerate(df_cols) if dc in norm_vars), None)
        if idx is not None: found[col_type] = df.columns[idx]
        else: return None
    return found

@st.cache_data(show_spinner=False)
def load_and_preprocess_data(uploaded_files, header_rows):
    st.info("Parsing files...")
    all_data = []
    for i, file in enumerate(uploaded_files):
        try:
            df = pd.read_excel(file, header=header_rows[i]-1, engine='openpyxl') if file.name.endswith('.xlsx') else pd.read_csv(file, header=header_rows[i]-1)
            df.columns = [c.strip() for c in df.columns]
            col_map = scan_for_columns(df)
            if col_map:
                df_clean = df.rename(columns={col_map['sample']: 'Sample Name', col_map['target']: 'Target Name', col_map['ct']: 'Ct'})
                df_clean['Ct'] = pd.to_numeric(df_clean['Ct'], errors='coerce')
                df_clean = df_clean[['Sample Name', 'Target Name', 'Ct']].dropna(subset=['Ct'])
                if not df_clean.empty:
                    all_data.append(df_clean)
                    st.success(f"✅ {file.name}")
        except Exception as e:
            st.error(f"Error: {file.name} - {e}")
    if not all_data: return None, None, None, None
    combined = pd.concat(all_data, ignore_index=True)
    combined['Sample Name'] = combined['Sample Name'].astype(str)
    combined['Target Name'] = combined['Target Name'].astype(str)
    combined['Group'] = combined['Sample Name'].str.strip()
    return combined, combined['Sample Name'].unique().tolist(), combined['Target Name'].unique().tolist(), combined['Group'].unique().tolist()

def calculate_ddct(data, control_gene, reference_group):
    ct_mean = data.groupby(['Group', 'Target Name'])['Ct'].mean().reset_index().rename(columns={'Ct': 'Mean Ct'})
    control_cts = ct_mean[ct_mean['Target Name'] == control_gene].set_index('Group')['Mean Ct']
    if control_cts.empty: return st.error(f"Control gene '{control_gene}' not found.")
    ct_mean['Control Ct'] = ct_mean['Group'].map(control_cts)
    ct_mean = ct_mean.dropna(subset=['Control Ct']).copy()
    ct_mean['ΔCt'] = ct_mean['Mean Ct'] - ct_mean['Control Ct']
    ref_dct = ct_mean[(ct_mean['Group'] == reference_group) & (ct_mean['Target Name'] != control_gene)].set_index('Target Name')['ΔCt']
    if ref_dct.empty: return st.error(f"Reference group '{reference_group}' invalid.")
    ct_mean['Reference ΔCt'] = ct_mean.apply(lambda r: ref_dct.get(r['Target Name']), axis=1)
    summary = ct_mean[ct_mean['Target Name'] != control_gene].dropna(subset=['Reference ΔCt']).copy()
    summary['ΔΔCt'] = summary['ΔCt'] - summary['Reference ΔCt']
    summary['Fold Change'] = 2**(-summary['ΔΔCt'])
    agg = data.groupby(['Group', 'Target Name'])['Ct'].agg(['std', 'count']).reset_index()
    control_agg = agg[agg['Target Name'] == control_gene].set_index('Group')[['std', 'count']].rename(columns={'std': 'Std Control', 'count': 'N Control'})
    full_agg = agg.merge(control_agg, on='Group', how='left')[agg['Target Name'] != control_gene].dropna(subset=['Std Control']).copy()
    full_agg['SE ΔCt'] = np.sqrt((full_agg['std']**2 / full_agg['count']) + (full_agg['Std Control']**2 / full_agg['N Control']))
    ref_se = full_agg[full_agg['Group'] == reference_group].set_index('Target Name')['SE ΔCt']
    full_agg['Ref SE'] = full_agg['Target Name'].map(ref_se)
    full_agg = full_agg.dropna(subset=['Ref SE']).copy()
    full_agg['SE ΔΔCt'] = np.sqrt(full_agg['SE ΔCt']**2 + full_agg['Ref SE']**2)
    summary = summary.merge(full_agg[['Group', 'Target Name', 'SE ΔΔCt']], on=['Group', 'Target Name'], how='left').set_index(['Group', 'Target Name'])
    summary['Upper'] = 2**(-(summary['ΔΔCt'] - summary['SE ΔΔCt']))
    summary['Lower'] = 2**(-(summary['ΔΔCt'] + summary['SE ΔΔCt']))
    p_vals = {}
    for (grp, tgt), _ in summary.iterrows():
        if grp == reference_group: p_vals[(grp, tgt)] = np.nan; continue
        try:
            g_vals = []
            for s in data[data['Group'] == grp]['Sample Name'].unique():
                tgt_ct = data[(data['Sample Name'] == s) & (data['Target Name'] == tgt)]['Ct'].values
                ctl_ct = data[(data['Sample Name'] == s) & (data['Target Name'] == control_gene)]['Ct'].values
                if len(tgt_ct) > 0 and len(ctl_ct) > 0:
                    val = float(tgt_ct[0]) - float(ctl_ct[0])
                    if not pd.isna(val): g_vals.append(val)
            r_vals = []
            for s in data[data['Group'] == reference_group]['Sample Name'].unique():
                tgt_ct = data[(data['Sample Name'] == s) & (data['Target Name'] == tgt)]['Ct'].values
                ctl_ct = data[(data['Sample Name'] == s) & (data['Target Name'] == control_gene)]['Ct'].values
                if len(tgt_ct) > 0 and len(ctl_ct) > 0:
                    val = float(tgt_ct[0]) - float(ctl_ct[0])
                    if not pd.isna(val): r_vals.append(val)
            p_vals[(grp, tgt)] = stats.ttest_ind(g_vals, r_vals, equal_var=False)[1] if len(g_vals) > 1 and len(r_vals) > 1 else np.nan
        except: p_vals[(grp, tgt)] = np.nan
    summary['P-Value'] = pd.Series(p_vals)
    summary[['Significant', 'Symbol']] = summary['P-Value'].apply(lambda p: pd.Series(('N/A', '') if pd.isna(p) else (('Yes (p<0.001)', '***') if p < 0.001 else (('Yes (p<0.01)', '**') if p < 0.01 else (('Yes (p<0.05)', '*') if p < 0.05 else ('No', 'ns'))))))
    return summary.reset_index()

def create_graph(summary, gene, cfg):
    plot_data = summary[summary['Target Name'] == gene].copy()
    if plot_data.empty: return None
    ref_grp = plot_data[plot_data['Fold Change'] == 1]['Group'].iloc[0] if not plot_data[plot_data['Fold Change'] == 1].empty else None
    plot_data['Error'] = plot_data['Upper'] - plot_data['Fold Change']
    plot_data.loc[plot_data['Group'] == ref_grp, 'Error'] = 0.0
    colors = [cfg['ref_color'] if r['Group'] == ref_grp else (cfg['sig_color'] if 'Yes' in r['Significant'] else cfg['nonsig_color']) for _, r in plot_data.iterrows()]
    fig, ax = plt.subplots(figsize=(cfg['width'], cfg['height']))
    ax.bar(plot_data['Group'], plot_data['Fold Change'], yerr=plot_data['Error'] if cfg['show_error_bars'] else None, capsize=cfg['capsize'], color=colors, edgecolor=cfg['edge_color'], linewidth=cfg['edge_width'], alpha=cfg['bar_alpha'])
    if cfg['show_significance']:
        for i, (_, r) in enumerate(plot_data.iterrows()):
            if 'Yes' in r['Significant']: ax.text(i, r['Upper'] + plot_data['Upper'].max() * 0.05, r['Symbol'], ha='center', va='bottom', fontsize=cfg['sig_fontsize'], color=cfg['sig_marker_color'], fontweight='bold')
    if cfg['show_ref_line']: ax.axhline(1, color=cfg['ref_line_color'], linestyle=cfg['ref_line_style'], linewidth=cfg['ref_line_width'], label='Reference')
    ax.set_ylabel(cfg['ylabel'], fontsize=cfg['label_fontsize'], fontweight='bold')
    ax.set_xlabel(cfg['xlabel'], fontsize=cfg['label_fontsize'], fontweight='bold')
    ax.set_title(cfg['title'].replace('{gene}', gene), fontsize=cfg['title_fontsize'], fontweight='bold')
    ax.set_ylim(0, plot_data['Upper'].max() * 1.25) if cfg['auto_ylim'] else ax.set_ylim(cfg['ymin'], cfg['ymax'])
    plt.xticks(rotation=cfg['xrotation'], ha='right')
    ax.spines['right'].set_visible(not cfg['remove_spines'])
    ax.spines['top'].set_visible(not cfg['remove_spines'])
    if cfg['show_grid']: ax.grid(axis='y', alpha=0.3, linestyle='--')
    if cfg['show_legend']: ax.legend(loc=cfg['legend_position'])
    plt.tight_layout()
    return fig

def create_excel(dfs_dict):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        for sheet, df in dfs_dict.items(): df.to_excel(writer, sheet_name=sheet, index=False)
    return output.getvalue()

def main():
    st.title("🧬 Interactive qPCR Analysis")
    st.sidebar.header("1. Upload Files")
    files = st.sidebar.file_uploader("Upload qPCR files", type=['csv', 'xlsx'], accept_multiple_files=True)
    if not files: return st.info("Upload files to begin.")
    st.sidebar.header("2. Header Rows")
    headers = [st.sidebar.number_input(f"{f.name}", 1, 100, 1, key=f'h{i}') for i, f in enumerate(files)]
    data, samples, targets, groups = load_and_preprocess_data(files, headers)
    if data is None: return
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    st.sidebar.header("📥 Exports")
    st.sidebar.download_button("📥 Raw Data", create_excel({'Raw': data}), f"raw_{ts}.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    
    st.header("📝 Step 3: Edit Samples")
    uniq = sorted([str(g) for g in groups])
    if 'sc' not in st.session_state: st.session_state['sc'] = {s: {'new': s, 'inc': True, 'ord': i} for i, s in enumerate(uniq)}
    t1, t2, t3 = st.tabs(["List", "Bulk", "Preview"])
    with t1:
        c1, c2, c3, c4 = st.columns([3,3,1,1])
        c1.markdown("**Original**"); c2.markdown("**New**"); c3.markdown("**Order**"); c4.markdown("**Include**")
        st.markdown("---")
        for s in uniq:
            c1, c2, c3, c4 = st.columns([3,3,1,1])
            cfg = st.session_state['sc'][s]
            c1.text(s)
            st.session_state['sc'][s]['new'] = c2.text_input("N", cfg['new'], key=f"n{s}", label_visibility="collapsed")
            st.session_state['sc'][s]['ord'] = c3.number_input("O", 0, len(uniq)-1, cfg['ord'], key=f"o{s}", label_visibility="collapsed")
            st.session_state['sc'][s]['inc'] = c4.checkbox("I", cfg['inc'], key=f"i{s}", label_visibility="collapsed")
    with t2:
        c1, c2 = st.columns(2)
        if c1.button("🔄 Reset"): st.session_state['sc'] = {s: {'new': s, 'inc': True, 'ord': i} for i, s in enumerate(uniq)}; st.rerun()
        if c1.button("✅ All"): [st.session_state['sc'][s].update({'inc': True}) for s in uniq]; st.rerun()
        if c1.button("❌ None"): [st.session_state['sc'][s].update({'inc': False}) for s in uniq]; st.rerun()
        pfx = c2.text_input("Prefix:", key="pfx")
        if c2.button("Add Prefix") and pfx: [st.session_state['sc'][s].update({'new': pfx + st.session_state['sc'][s]['new']}) for s in uniq]; st.rerun()
        sfx = c2.text_input("Suffix:", key="sfx")
        if c2.button("Add Suffix") and sfx: [st.session_state['sc'][s].update({'new': st.session_state['sc'][s]['new'] + sfx}) for s in uniq]; st.rerun()
    with t3:
        sorted_s = sorted(st.session_state['sc'].items(), key=lambda x: x[1]['ord'])
        prev = [{'Order': c['ord'], 'Original': o, 'New': c['new'], 'Status': '✅' if c['inc'] else '❌'} for o, c in sorted_s]
        c1, c2, c3 = st.columns(3)
        c1.metric("Total", len(uniq))
        c2.metric("Included", sum(1 for c in st.session_state['sc'].values() if c['inc']))
        c3.metric("Excluded", len(uniq) - sum(1 for c in st.session_state['sc'].values() if c['inc']))
        st.dataframe(pd.DataFrame(prev), use_container_width=True, hide_index=True)
    
    if st.button("✅ Apply Changes", type="primary"):
        inc_orig = [o for o, c in st.session_state['sc'].items() if c['inc']]
        if not inc_orig: return st.error("Select at least one sample!")
        filt = data[data['Group'].isin(inc_orig)].copy()
        name_map = {o: c['new'] for o, c in st.session_state['sc'].items()}
        filt['Group'] = filt['Group'].map(name_map)
        filt['Sample Name'] = filt['Sample Name'].map(lambda x: name_map.get(x, x))
        grps_filt = [c['new'] for o, c in sorted(st.session_state['sc'].items(), key=lambda x: x[1]['ord']) if c['inc']]
        st.session_state.update({'df': filt, 'gf': grps_filt, 'ec': True})
        st.success("✅ Applied!"); st.rerun()
    
    if 'ec' not in st.session_state or not st.session_state['ec']: return st.info("👆 Apply changes first.")
    data, groups = st.session_state['df'], st.session_state['gf']
    st.sidebar.download_button("📥 Edited Data", create_excel({'Edited': data}), f"edited_{ts}.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    
    st.sidebar.header("4. Analysis")
    cg = st.sidebar.selectbox("Control Gene:", targets, index=targets.index('ACTIN') if 'ACTIN' in targets else 0)
    rg = st.sidebar.selectbox("Reference Group:", groups, index=0)
    if st.sidebar.button("🔬 Calculate"):
        summ = calculate_ddct(data, cg, rg)
        if summ is not None: st.session_state.update({'summ': summ, 'cg': cg}); st.success("✅ Done!")
    
    if 'summ' not in st.session_state: return
    summ, cg = st.session_state['summ'], st.session_state['cg']
    st.header("📊 Results")
    tgts = [t for t in summ['Target Name'].unique() if t != cg]
    for g in tgts:
        with st.expander(f"🧬 {g}", expanded=True):
            gd = summ[summ['Target Name'] == g].copy()
            for c in ['Mean Ct', 'ΔCt', 'ΔΔCt', 'Fold Change', 'SE ΔΔCt', 'Upper', 'Lower', 'P-Value']:
                if c in gd.columns: gd[c] = gd[c].round(3)
            st.dataframe(gd, hide_index=True, use_container_width=True)
    st.sidebar.download_button("📥 Results", create_excel({'Summary': summ}), f"results_{ts}.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    
    st.header("📈 Graphs")
    with st.expander("🎨 Customize", expanded=True):
        c1, c2, c3 = st.columns(3)
        w, h, xr = c1.slider("Width", 6, 16, 10), c1.slider("Height", 4, 12, 6), c1.slider("X rotation", 0, 90, 45)
        rc, sc, nc, ec = c2.color_picker("Ref", "#95a5a6"), c2.color_picker("Sig", "#e74c3c"), c2.color_picker("NonSig", "#3498db"), c2.color_picker("Edge", "#000000")
        seb, ss, srl, sg, sl = c3.checkbox("Error bars", True), c3.checkbox("Sig markers", True), c3.checkbox("Ref line", True), c3.checkbox("Grid", True), c3.checkbox("Legend", True)
        c1, c2, c3 = st.columns(3)
        ttl, xl, yl = c1.text_input("Title", "Expression: {gene}"), c1.text_input("X-label", "Group"), c1.text_input("Y-label", "Fold Change")
        tfs, lfs, sfs = c2.slider("Title size", 10, 20, 14), c2.slider("Label size", 8, 16, 12), c2.slider("Sig size", 12, 24, 18)
        ba, ew, cs = c3.slider("Alpha", 0.0, 1.0, 0.8), c3.slider("Edge width", 0.0, 3.0, 1.2), c3.slider("Capsize", 0, 10, 5)
        c1, c2 = st.columns(2)
        ay = c1.checkbox("Auto Y", True)
        ymn, ymx = (c1.number_input("Y min", value=0.0), c1.number_input("Y max", value=3.0)) if not ay else (0, 3)
        rlc, rls, rlw = c1.color_picker("Ref line", "#808080"), c1.selectbox("Line style", ['--', '-', '-.', ':']), c1.slider("Ref width", 0.5, 3.0, 1.0)
        smc, rs, lp = c2.color_picker("Sig color", "#e74c3c"), c2.checkbox("Remove spines", True), c2.selectbox("Legend pos", ['upper left', 'upper right', 'lower left', 'lower right', 'best'], 1)
    
    cfg = {'width': w, 'height': h, 'xrotation': xr, 'ref_color': rc, 'sig_color': sc, 'nonsig_color': nc, 'edge_color': ec, 'show_error_bars': seb, 'show_significance': ss, 'show_ref_line': srl, 'show_grid': sg, 'show_legend': sl, 'title': ttl, 'xlabel': xl, 'ylabel': yl, 'title_fontsize': tfs, 'label_fontsize': lfs, 'sig_fontsize': sfs, 'bar_alpha': ba, 'edge_width': ew, 'capsize': cs, 'auto_ylim': ay, 'ymin': ymn, 'ymax': ymx, 'ref_line_color': rlc, 'ref_line_style': rls, 'ref_line_width': rlw, 'sig_marker_color': smc, 'remove_spines': rs, 'legend_position': lp}
    sel = st.multiselect("Select genes:", tgts, default=tgts)
    if st.button("🎨 Generate", type="primary"):
        if sel:
            tabs = st.tabs([f"📊 {g}" for g in sel])
            zb = io.BytesIO()
            with zipfile.ZipFile(zb, 'w', zipfile.ZIP_DEFLATED) as zf:
                for i, g in enumerate(sel):
                    with tabs[i]:
                        fig = create_graph(summ, g, cfg)
                        if fig:
                            st.pyplot(fig)
                            ib = io.BytesIO()
                            fig.savefig(ib, format='png', dpi=300, bbox_inches='tight')
                            ib.seek(0)
                            st.download_button(f"📥 {g}", ib.getvalue(), f"{g}_{ts}.png", "image/png")
                            ib.seek(0)
                            zf.writestr(f"{g}.png", ib.getvalue())
                            plt.close(fig)
            zb.seek(0)
            st.download_button("📥 All Graphs (ZIP)", zb.getvalue(), f"graphs_{ts}.zip", "application/zip")
        else: st.warning("Select at least one gene.")
    
    st.header("📋 Statistics")
    for g in tgts:
        with st.expander(f"📊 {g}"):
            gd = summ[summ['Target Name'] == g]
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Mean FC", f"{gd['Fold Change'].mean():.3f}")
            c2.metric("Significant", f"{len(gd[gd['Significant'].str.contains('Yes')])}/{len(gd)}")
            c3.metric("Max FC", f"{gd['Fold Change'].max():.3f}")
            c4.metric("Min P", f"{gd['P-Value'].min():.4f}" if pd.notna(gd['P-Value'].min()) else "N/A")

if __name__ == "__main__": main()
