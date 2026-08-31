import re, unicodedata
from pathlib import Path
import numpy as np, pandas as pd, streamlit as st, matplotlib.pyplot as plt
import streamlit.components.v1 as components

st.set_page_config(page_title='Protein-oriented HOF prioritization',layout='wide',initial_sidebar_state='expanded')
st.markdown('''<style>:root{color-scheme:light!important}html,body,.stApp{background:#fff!important;color:#111!important}section[data-testid="stSidebar"]{background:#f6f7f8!important}section[data-testid="stSidebar"] *{color:#111!important}div[data-baseweb="select"] *,div[role="listbox"] *{color:#111!important;-webkit-text-fill-color:#111!important}</style>''',unsafe_allow_html=True)
P='protein_descriptors_FINAL_100_UI_MIN.csv'; W='HOF_MASTER_WORKBOOK_functional_group_typed_FIXED.xlsx'; C='HOF_DATABASE_FINAL_CURATED_V1.xlsx'; D=Path('hof_db')
for f in [P,W,C]:
    if not Path(f).exists(): st.error(f'Missing required file: {f}'); st.stop()
def n(x): return pd.to_numeric(x,errors='coerce')
def mm(x):
    s=n(x); a,b=s.min(),s.max(); return pd.Series(0.,index=s.index) if pd.isna(a) or pd.isna(b) or b-a<=1e-15 else ((s-a)/(b-a)).fillna(0.)
def pc(x): return n(x).rank(method='average',pct=True).fillna(0.)
def fk(x): return re.sub(r'[^a-z0-9]+','',unicodedata.normalize('NFKC',str(x)).casefold().replace('.cif',''))
def fm(x,k=3):
    try: return 'NA' if pd.isna(x) else f'{float(x):.{k}f}'
    except: return 'NA'
def fclass(r):
    fn=str(r.get('cif_filename','')); c=int(r.ccdc_id) if pd.notna(r.get('ccdc_id')) else None; m=float(r.get('metal_count_cif',0) or 0) if pd.notna(r.get('metal_count_cif')) else 0.; t=(fn+' '+str(r.get('discovery_publication_citation') or '')).lower()
    if c==1991981:return 'retained_verified_hof'
    if c in {1991980,1991982,1566036,1566037,1893340,1893341}:return 'excluded_confirmed_nonhof'
    if fn=='257-rub-hofs.cif':return 'retained_verified_metal_hof'
    if m<=0:return 'retained_source_hof_candidate'
    if any(z in t for z in ['hydrogen bonded','hydrogen-bonded','m-hof','hof-','hofs']):return 'retained_verified_metal_hof'
    if any(z in t for z in ['coordination polymer','metal-organic framework','porous 3d mof',' mofs ']):return 'excluded_confirmed_nonhof'
    return 'manual_review_excluded_main'

@st.cache_data(show_spinner='Building audited 651 × 100 matrix…')
def load():
    p=pd.read_csv(P); h=pd.read_excel(W,sheet_name='core_with_functional_groups').drop_duplicates('cif_filename'); s=pd.read_excel(C,sheet_name='framework_series_layer')
    assert len(p)==100 and p.uniprot_id.nunique()==100
    p['protein_id_eval']=p.uniprot_id.astype(str)
    for x in ['surface_hbond_donor_density','surface_hbond_acceptor_density','surface_pos_residue_frac','surface_neg_residue_frac','surface_hydrophobic_frac','surface_aromatic_frac']: p[x+'_norm']=mm(p[x])
    h['framework_class_status']=h.apply(fclass,axis=1); h['main_analysis_included']=h.framework_class_status.str.startswith('retained')
    sk=[x for x in ['cif_filename','framework_series_key','framework_series','family_revised','broad_family_final','display_group_final','series_size','representative_cif','is_representative'] if x in s.columns]
    h=h.merge(s[sk].drop_duplicates('cif_filename'),on='cif_filename',how='left'); miss=h.framework_series_key.isna()|h.framework_series_key.astype(str).str.strip().eq('')
    h['series_key_eval']=h.framework_series_key; h.loc[miss,'series_key_eval']=h.loc[miss,'cif_filename'].map(fk); h['series_label_eval']=h.framework_series.where(h.framework_series.notna(),h.series_key_eval)
    q=h.ccdc_id.eq(1991981); h.loc[q,'series_key_eval']='pentahof1'; h.loc[q,'series_label_eval']='pentaHOF-1'
    for x in ['Di_A_best','Df_A_best','Dif_A_best','ASA_m2_g_best','AV_volume_fraction_best']: h[x]=n(h[x]).fillna(0.)
    di=h.Di_A_best; h['relative_aperture_path_ratio']=.5*np.where(di>0,(h.Df_A_best/di).clip(0,1),0)+.5*np.where(di>0,(h.Dif_A_best/di).clip(0,1),0)
    h=h[h.main_analysis_included & h.fg_typing_success.eq(True)].copy(); assert len(h)==651
    for x in ['hbond_donor_group_count_per_100_heavy_atoms_scaled','hbond_acceptor_group_count_per_100_heavy_atoms_scaled','raw_charge_motif_score_scaled','raw_aromatic_motif_score_scaled']: h[x]=n(h[x]).fillna(0).clip(0,1)
    h['Di_norm']=mm(h.Di_A_best);h['AV_norm']=mm(h.AV_volume_fraction_best);h['ASA_norm']=mm(h.ASA_m2_g_best);h['O_norm']=mm(h.relative_aperture_path_ratio)
    h['G_enc']=.40*h.Di_norm+.35*h.AV_norm+.25*h.O_norm; h['G_surf']=.70*h.ASA_norm+.30*h.O_norm; h['P_G_enc']=pc(h.G_enc);h['P_G_surf']=pc(h.G_surf)
    hc=['cif_filename','series_key_eval','series_label_eval','Di_A_best','Df_A_best','Dif_A_best','ASA_m2_g_best','AV_volume_fraction_best','G_enc','G_surf','P_G_enc','P_G_surf','protein_screen_confidence_tier','hbond_donor_group_count_per_100_heavy_atoms_scaled','hbond_acceptor_group_count_per_100_heavy_atoms_scaled','raw_charge_motif_score_scaled','raw_aromatic_motif_score_scaled']
    pr=['protein_id_eval','uniprot_id','protein_name','organism','min_dimension_A','effective_diameter_A','mean_plddt','surface_residue_fraction','surface_hbond_donor_density_norm','surface_hbond_acceptor_density_norm','surface_pos_residue_frac_norm','surface_neg_residue_frac_norm','surface_hydrophobic_frac_norm','surface_aromatic_frac_norm']
    a=h[hc].copy();b=p[pr].copy();a['_k']=b['_k']=1;r=a.merge(b,on='_k').drop(columns='_k')
    r['H']=.5*(r.hbond_donor_group_count_per_100_heavy_atoms_scaled*r.surface_hbond_acceptor_density_norm+r.hbond_acceptor_group_count_per_100_heavy_atoms_scaled*r.surface_hbond_donor_density_norm)
    r['Q']=r.raw_charge_motif_score_scaled*.5*(r.surface_pos_residue_frac_norm+r.surface_neg_residue_frac_norm); r['R']=r.raw_aromatic_motif_score_scaled*.5*(r.surface_hydrophobic_frac_norm+r.surface_aromatic_frac_norm)
    r['chemistry']=.35*r.H+.35*r.Q+.30*r.R; r['Pp_C']=r.groupby('protein_id_eval').chemistry.rank(method='average',pct=True); r['S_hp']=r.groupby('cif_filename').chemistry.rank(method='average',pct=True)
    meta=[x for x in ['cif_filename','family','family_original','discovery_publication_year','discovery_publication_doi','discovery_publication_citation','final_doi','fg_confidence','fg_flags','family_revised','broad_family_final','display_group_final'] if x in h.columns]
    r=r.merge(h[meta].drop_duplicates('cif_filename'),on='cif_filename',how='left'); assert len(r)==65100
    return r,p
R,PANEL=load()

@st.cache_data(show_spinner=False)
def vindex():
    z={}
    if D.exists():
        for d in D.iterdir():
            if d.is_dir(): z.setdefault(fk(d.name),[]).append(d)
    return z
VI=vindex()
def views(c):
    ds=VI.get(fk(c),[])
    if len(ds)!=1:return {}
    return {x:str(ds[0]/f'{x}.html') for x in ['unit','super','wire'] if (ds[0]/f'{x}.html').exists()}
def lit(r):
    return [(x,str(r.get(x))) for x in ['final_doi','discovery_publication_doi','discovery_publication_citation'] if pd.notna(r.get(x)) and str(r.get(x)).strip()]
def family(r):
    for x in ['broad_family_final','family_revised','family','family_original','display_group_final']:
        if pd.notna(r.get(x)) and str(r.get(x)).strip():return str(r.get(x))
    return 'Unclassified'
def route(pid,kind,t):
    d=R[R.protein_id_eval.eq(pid)].copy()
    if kind=='Growth-mediated integration / encapsulation': d['PG']=d.P_G_enc; g=(d.AV_volume_fraction_best>0)&(d.ASA_m2_g_best>0)&(d.P_G_enc>=t);note='Framework forms around/with the protein; final pore aperture is not treated as a whole-protein entry requirement.'
    elif kind=='Accessible-interface contact': d['PG']=d.P_G_surf;g=(d.ASA_m2_g_best>0)&(d.P_G_surf>=t);note='Small-probe-accessible framework interfaces are compared; periodic ASA is not claimed to equal external particle area.'
    else: d['PG']=d.P_G_enc;g=(d.AV_volume_fraction_best>0)&(d.ASA_m2_g_best>0)&(d.P_G_enc>=t)&(d.Df_A_best>=d.min_dimension_A)&(d.Di_A_best>=d.effective_diameter_A);note='Conservative static whole-protein entry filter; short or empty lists are expected.'
    d=d[g].copy();d['final_score']=.65*d.Pp_C+.30*d.S_hp+.05*d.PG;return d,note

def profile(r):
    labs=['H','Q','R','Pp_C','S_hp','Geometry'];v=[r.H,r.Q,r.R,r.Pp_C,r.S_hp,r.PG];a=np.linspace(0,2*np.pi,len(labs),endpoint=False);a=np.r_[a,a[0]];v=np.r_[v,v[0]];fig=plt.figure(figsize=(4,4));ax=fig.add_subplot(111,polar=True);ax.plot(a,v);ax.fill(a,v,alpha=.12);ax.set_xticks(a[:-1]);ax.set_xticklabels(labs);ax.set_ylim(0,1);return fig

st.sidebar.header('Screening controls'); opts=['— Select protein —']+[f'{x.protein_name} | {x.uniprot_id}' for _,x in PANEL.sort_values(['protein_name','uniprot_id']).iterrows()];sel=st.sidebar.selectbox('Target protein',opts)
kind=st.sidebar.selectbox('Screening route',['Growth-mediated integration / encapsulation','Accessible-interface contact','Strict post-synthetic infiltration']);tau=st.sidebar.slider('Geometry feasibility threshold (τ)',0.,.75,.25,.05);top=st.sidebar.slider('Recommendations',5,50,20,5);collapse=st.sidebar.checkbox('Collapse framework-series variants',True);r3=st.sidebar.checkbox('Require exact 3D view');rl=st.sidebar.checkbox('Require literature metadata');q=st.sidebar.text_input('Search CIF / series / family','')
st.sidebar.markdown('---');st.sidebar.caption('Frozen manuscript model');st.sidebar.code('H/Q/R = 0.35/0.35/0.30\nFinal = 0.65 Pp_C + 0.30 S_hp + 0.05 P_G');st.sidebar.caption('Literature, family labels and 3D availability never enter the score.')
st.title('Protein-oriented HOF prioritization');st.markdown('**Auditable decision support for the frozen 100-protein / 651-HOF analysis.** Scores are prioritization hypotheses, not experimental compatibility probabilities.')
a,b,c,d=st.columns(4);a.metric('Curated proteins',100);b.metric('Typed HOFs',651);c.metric('Pairwise chemistry rows','65,100');d.metric('Typed framework series',550)
if sel=='— Select protein —':st.info('Select a protein to generate a route-specific shortlist.');st.stop()
pr=PANEL[(PANEL.protein_name+' | '+PANEL.uniprot_id).eq(sel)].iloc[0];pid=str(pr.uniprot_id);df,note=route(pid,kind,tau)
if r3:df=df[df.cif_filename.map(lambda x:bool(views(x)))]
if rl:df=df[df.apply(lambda x:bool(lit(x)),axis=1)]
if q.strip():
    z=q.casefold();df=df[df.cif_filename.astype(str).str.casefold().str.contains(z,na=False)|df.series_label_eval.astype(str).str.casefold().str.contains(z,na=False)|df.apply(lambda x:z in family(x).casefold(),axis=1)]
df=df.sort_values(['final_score','PG','cif_filename'],ascending=[False,False,True]);
if collapse:df=df.drop_duplicates('series_key_eval')
df=df.head(top).reset_index(drop=True);df.insert(0,'rank',np.arange(1,len(df)+1))
st.subheader(f'{pr.protein_name} — {pid}');st.caption(note);x1,x2,x3,x4=st.columns(4);x1.metric('Organism',str(pr.get('organism','NA'))[:28]);x2.metric('Mean pLDDT',fm(pr.get('mean_plddt'),1));x3.metric('Surface residues',int(pr.get('surface_residue_count',0)));x4.metric('Candidates shown',len(df))
if (pd.notna(pr.get('mean_plddt')) and float(pr.mean_plddt)<70) or (pd.notna(pr.get('surface_residue_fraction')) and float(pr.surface_residue_fraction)>=.95):st.warning('This protein belongs to the lower-confidence structural sensitivity subset; interpret surface-specific ranking cautiously.')
show=['rank','series_label_eval','cif_filename','final_score','Pp_C','S_hp','PG','H','Q','R','Df_A_best','Di_A_best','AV_volume_fraction_best','ASA_m2_g_best'];st.dataframe(df[show],use_container_width=True,hide_index=True);st.download_button('Download shortlist',df[show].to_csv(index=False),f'HOF_{pid}_top{top}.csv','text/csv')
for _,r in df.iterrows():
    st.markdown('---');st.markdown(f"### #{int(r['rank'])} — {r.series_label_eval}");st.caption(f"{r.cif_filename} · {family(r)} · HOF confidence: {r.get('protein_screen_confidence_tier','NA')}")
    y1,y2,y3,y4=st.columns(4);y1.metric('Final',fm(r.final_score));y2.metric('Chemistry %ile',fm(r.Pp_C));y3.metric('Differentiation',fm(r.S_hp));y4.metric('Geometry',fm(r.PG))
    t1,t2,t3,t4=st.tabs(['Why shortlisted','Score profile','3D view','Literature'])
    with t1:
        pros=[];cau=[]
        pros.append('Protein-specific chemistry is in the top quartile.' if r.Pp_C>=.75 else 'Candidate passes the selected route gate.');pros.append('Strong target-relative differentiation.' if r.S_hp>=.75 else 'More generalist than target-specific.' if r.S_hp<.5 else 'Moderate target-relative differentiation.');pros.append('Largest raw chemistry contribution: '+max({'H-bond':r.H,'polar/charged motif':r.Q,'aromatic/hydrophobic':r.R},key={'H-bond':r.H,'polar/charged motif':r.Q,'aromatic/hydrophobic':r.R}.get));
        if r.PG<.35:cau.append('Route geometry is relatively weak despite passing the threshold.')
        if kind=='Strict post-synthetic infiltration':cau.append('Static dimensions do not model flexibility, solvent or kinetics.')
        st.markdown('**Supporting evidence**');[st.markdown('- '+z) for z in pros];st.markdown('**Cautions**');[st.markdown('- '+z) for z in cau] if cau else st.markdown('- No additional descriptor-level caution triggered.');st.caption('Descriptor explanation only; not an experimental outcome prediction.')
    with t2:st.pyplot(profile(r),use_container_width=False)
    with t3:
        v=views(r.cif_filename)
        if not v:st.info('No unique exact-normalized 3D HTML view for this CIF.')
        else:
            z=st.radio('View',list(v),horizontal=True,key=f"v_{pid}_{r['rank']}");components.html(Path(v[z]).read_text(encoding='utf-8',errors='ignore'),height=620,scrolling=True)
    with t4:
        L=lit(r)
        if not L:st.info('No DOI/citation metadata attached to this structural entry.')
        else:
            for k,v in L:st.markdown(f'**{k}:** {v}')
        st.caption('Literature is traceability/context only and is not a score input.')
st.markdown('---');st.caption('Frozen analysis: 100 unique proteins, 651 functionally typed HOFs, 65,100 pairwise records; default τ = 0.25.')
