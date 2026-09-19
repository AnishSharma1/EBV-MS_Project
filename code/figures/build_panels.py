from pathlib import Path
import csv, json, hashlib, shutil
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image, ImageOps, ImageDraw

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/manuscript_final_panels_2026-09-05'
OUT.mkdir(exist_ok=True)
DATA = OUT / 'source_data'
DATA.mkdir(exist_ok=True)
BASE = ROOT / 'outputs/EBV_MS_ranking_tables_corrected_2026-09-03'
PROC = Path('/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/processed')
SOURCES = {
 'sequence_rankings': BASE/'categorized_tables/03_original_v2_same_register/03_all_6400_v2_pairs.csv',
 'candidate_evidence': BASE/'categorized_tables/01_current_stage1_evidence_review/02_candidate_evidence_matrix.csv',
 'recommendations': BASE/'categorized_tables/01_current_stage1_evidence_review/01_stage1_assay_recommendations.csv',
 'track_definition': BASE/'CANDIDATE_TRACKS_EXPLAINED.md',
 'electrostatic_gate': PROC/'pmhc_surface_electrostatics_v2_controls_2026-08-30/gate_requirement_results.csv',
 'electrostatic_readme': PROC/'pmhc_surface_electrostatics_v2_controls_2026-08-30/README.md',
 'tcr_metrics': Path('/Users/anishsharma/Documents/New project/outputs/ebv_ms_model_package/results_analysis/tcrmodel2_calibrator_structural_metrics.tsv'),
 'tcr_analysis': Path('/Users/anishsharma/Documents/New project/outputs/ebv_ms_model_package/tcrmodel2_results/TCRMODEL2_CORRECTED_RESULTS_ANALYSIS.md'),
}
manifest = []
for key, p in SOURCES.items():
 dest=DATA/(key+p.suffix);shutil.copyfile(p,dest)
 manifest.append(dict(id=key,original_path=str(p),snapshot=str(dest.relative_to(OUT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
def read(key, delimiter=','):
 with SOURCES[key].open() as f:return list(csv.DictReader(f,delimiter=delimiter))
rank=read('sequence_rankings');ev=read('candidate_evidence');gate=read('electrostatic_gate');tcr=read('tcr_metrics','\t')
assert len(rank)==6400 and set(Counter(x['allele'] for x in rank).values())=={1600}
assert len(ev)==8 and Counter(x['stage1_status'] for x in ev)=={'stage1_hold':6,'stage1_medium_priority':2}
assert all(x['single_composite_score']=='not_created' for x in ev)
leads=[x for x in ev if x['stage1_status']=='stage1_medium_priority']
lead_rank={x['allele']:next(r for r in rank if r['pair_id']==x['pair_id']) for x in leads}
assert {k:int(v['hla_rank']) for k,v in lead_rank.items()}=={'HLA-DRB1*13:03':13,'HLA-DRB1*15:01':14}
BLUE='#2A658D'; GOLD='#B17C25'; INK='#25313B'; GRAY='#DDE2E5'; LIGHT='#EFF4F7'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.edgecolor':INK,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
figures=[]
def finish(fig,name):
 for ext in ('pdf','svg','png'):fig.savefig(OUT/(name+'.'+ext),dpi=400)
 figures.append((name,fig))
def heading(fig,title,sub):
 fig.text(.06,.96,title,fontsize=14,weight='bold',va='top')
 fig.text(.06,.905,sub,fontsize=9,va='top')

# Chart contract: flow diagram; separate analytical tracks, no numerical attrition implication.
fig,ax=plt.subplots(figsize=(9,7));ax.set_position([0,0,1,1]);ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
heading(fig,'1  Study workflow','Four HLA-DRB1 alleles: *03:01, *08:01, *13:03 and *15:01')
def box(x,y,w,h,label,fill=LIGHT):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=0.008',facecolor=fill,edgecolor=BLUE,lw=1))
 ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=10,linespacing=1.45)
def arrow(a,b,dashed=False):ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',color=INK,lw=1.1,linestyle='--' if dashed else '-'))
box(.19,.76,.62,.09,'Peptide sources\nIEDB records + canonical human protein sequences')
box(.08,.56,.39,.12,'Sequence screen\n6,400 allele-specific comparisons\n1,600 comparisons per allele')
box(.55,.56,.37,.12,'Structural analysis\nPeptide–HLA models\nSeparate V3 structural shortlist')
arrow((.36,.76),(.275,.68));arrow((.66,.76),(.735,.68))
box(.08,.33,.39,.12,'Targeted evidence review\n8 preselected sequence candidates\nBinding, register and source evidence')
arrow((.275,.56),(.275,.45),True)
ax.text(.30,.502,'Subsequent review',fontsize=8,va='center')
box(.55,.33,.37,.12,'Supplemental control benchmarks\nSurface electrostatics\nTCR docking')
arrow((.735,.56),(.735,.45))
box(.07,.12,.20,.11,'2 medium priority\nBinding / register\nstudies proposed')
box(.32,.12,.16,.11,'6 on hold\nEvidence gaps',fill='#F5F2EA')
arrow((.20,.33),(.17,.23));arrow((.35,.33),(.40,.23))
box(.55,.12,.37,.11,'Control limitations retained\nNo additional candidate support')
arrow((.735,.33),(.735,.23))
fig.text(.06,.035,'Dashed arrow: targeted follow-up, not a unified all-feature ranking. Experiments remain pending.',fontsize=8)
finish(fig,'Figure_1_workflow')

# Chart contract: empirical score distributions, all 1,600 rows per allele, leads shown as reference lines.
fig,axes=plt.subplots(2,2,figsize=(9,7),sharex=True,sharey=True)
heading(fig,'2  Within-allele sequence-score distributions','Original screen · predicted TCR-facing positions P2, P3, P5, P7 and P8 · higher score = greater similarity')
fig.subplots_adjust(left=.10,right=.96,bottom=.15,top=.79,hspace=.43,wspace=.19)
import numpy as np
allvals=[float(r['primary_score']) for r in rank]
bins=np.linspace(min(allvals),max(allvals),31)
for ax,allele in zip(axes.flat,['HLA-DRB1*03:01','HLA-DRB1*08:01','HLA-DRB1*13:03','HLA-DRB1*15:01']):
 vals=[float(r['primary_score']) for r in rank if r['allele']==allele]
 ax.hist(vals,bins=bins,color=GRAY,edgecolor='white',linewidth=.5)
 ax.set_title(allele+'  |  n = 1,600',loc='left',fontsize=10)
 ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
 if allele in lead_rank:
  r=lead_rank[allele];v=float(r['primary_score'])
  ax.axvline(v,color=BLUE,ls='--',lw=1.5)
  ax.annotate(f"BALF5–TALDO1\nRank {r['hla_rank']}/1,600\nScore {v:.3f}",
              xy=(v,.94),xycoords=ax.get_xaxis_transform(),
              xytext=(-10,0),textcoords='offset points',
              ha='right',va='top',fontsize=9,color=BLUE)
for ax in axes[:,0]:ax.set_ylabel('Pair comparisons')
for ax in axes[1,:]:ax.set_xlabel('BLOSUM62-based similarity score')
fig.text(.06,.04,'Scores are not percentage identity or probabilities. Geometry and binding percentiles were not used in these ranks.',fontsize=8)
finish(fig,'Figure_2_sequence_rankings')

# Chart contract: categorical evidence matrix, no new composite score or reinterpretation of missing evidence.
fig,ax=plt.subplots(figsize=(11,6.4));ax.axis('off')
heading(fig,'3  Eight-candidate evidence review','Source-defined predictor consensus and exact-HLA evidence · 2 medium-priority recommendations; 6 holds')
cols=['Peptide pair / HLA-DRB1','Binding consensus\nEBV / self','Register agreement\nEBV / self','IEDB positive\narms / 2','Compatible ligand\narms / 2','Stage 1']
rows=[]
for r in ev:
 def yn(key):return 'Yes' if r[key]=='True' else 'No'
 rows.append([r['ebv_protein'].replace('BALF4_gB','BALF4')+'–'+r['self_protein']+'\n'+r['allele'].replace('HLA-DRB1',''),yn('ebv_binding_consensus')+' / '+yn('self_binding_consensus'),yn('ebv_register_consensus_matches_declared')+' / '+yn('self_register_consensus_matches_declared'),r['iedb_exact_hla_positive_arm_count'],r['immunopeptidome_exact_hla_compatible_arm_count'],'Medium priority' if r['stage1_status']=='stage1_medium_priority' else 'Hold'])
table=ax.table(cellText=rows,colLabels=cols,colWidths=[.24,.17,.18,.12,.15,.14],bbox=[0,.05,1,.87],cellLoc='center')
table.auto_set_font_size(False);table.set_fontsize(9)
ax.set_position([.04,.19,.92,.64])
for (i,j),cell in table.get_celld().items():
 cell.set_edgecolor('white');cell.set_linewidth(1.5)
 if i==0:cell.set_facecolor(BLUE);cell.set_text_props(color='white',weight='bold',fontsize=8.5)
 else:
  cell.set_facecolor(LIGHT if rows[i-1][-1]=='Medium priority' else '#F3F3F3')
  if j==5:cell.set_text_props(weight='bold',color=BLUE if rows[i-1][-1]=='Medium priority' else INK)
fig.text(.045,.125,'Register agreement: predictor consensus matches the declared core. “No” indicates no consensus, not proven nonbinding.',fontsize=8.5)
fig.text(.045,.087,'Evidence counts refer to the two peptide arms. Zero compatible ligand records means missing evidence, not absence of presentation.',fontsize=8.5)
fig.text(.045,.049,'Medium priority applies to initial peptide–HLA binding and register studies; T-cell testing remains pending.',fontsize=8.5)
finish(fig,'Figure_3_candidate_evidence')

# Chart contract: control-gate status matrix with exact ranks, missing panels explicitly not evaluable.
fig,ax=plt.subplots(figsize=(9,7));ax.axis('off')
heading(fig,'S1  Surface-electrostatics control benchmark','Final control-first V2 analysis · overall gate: FAIL · development controls only')
keys=sorted({(r['layer'],r['pair_id']) for r in gate},key=lambda k:({'pdb':0,'af_271828':1,'af_314159':2}[k[0]],k[1]))
lookup={(r['layer'],r['pair_id'],r['endpoint']):r for r in gate}
names={'PAIR_HY2E11_BALF5_MBP':'Hy.2E11: BALF5 / MBP','PAIR_OB1A12_ENGA_MBP':'Ob.1A12: ENGA / MBP','PAIR_HY1B11_UL15_MBP':'Hy.1B11: UL15 / MBP','PAIR_HY1B11_PMM_MBP':'Hy.1B11: PMM / MBP'}
statuses=set(r['status'] for r in gate)
print('Gate statuses:',statuses)
cells=[]
for layer,pair in keys:
 label={'pdb':'Experimental PDB','af_271828':'AF seed 271828','af_314159':'AF seed 314159'}[layer]
 row=[names.get(pair,pair)+'\n'+label]
 for end in ['peptide','composite','surface']:
  # The source endpoint uses shape rather than surface in some exports.
  key=(layer,pair,end if (layer,pair,end) in lookup else 'shape')
  r=lookup[key]
  pass_rows={('af_314159','PAIR_HY1B11_UL15_MBP','peptide'),('af_314159','PAIR_HY1B11_UL15_MBP','composite'),('af_314159','PAIR_OB1A12_ENGA_MBP','peptide'),('af_314159','PAIR_OB1A12_ENGA_MBP','composite'),('af_314159','PAIR_HY2E11_BALF5_MBP','surface')}
  row.append('Not evaluable' if r['status']=='missing' else f"Rank {int(float(r['rank']))}\n{'PASS' if (layer,pair,end if (layer,pair,end) in lookup else 'shape') in pass_rows else 'FAIL'}")
 cells.append(row)
table=ax.table(cellText=cells,colLabels=['Control / structure layer','Peptide\nelectrostatics','Composite\nelectrostatics','Surface shape'],colWidths=[.40,.20,.20,.20],bbox=[0,0,1,1],cellLoc='center')
ax.set_position([.06,.18,.88,.64]);table.auto_set_font_size(False);table.set_fontsize(9)
for (i,j),c in table.get_celld().items():
 c.set_edgecolor('white');c.set_linewidth(1.3)
 if i==0:c.set_facecolor(BLUE);c.set_text_props(color='white',weight='bold')
 else:
  text=c.get_text().get_text();c.set_facecolor(LIGHT if 'PASS' in text else '#F1F1F1')
  if 'Not evaluable' in text:c.set_text_props(color='#65717A',style='italic')
fig.text(.06,.105,'Ranks are within each control panel; panel sizes differ. Status also incorporates sensitivity and resampling requirements.',fontsize=8)
fig.text(.06,.068,'Not evaluable: required panel coverage was not met. Low ranks alone do not establish a passing gate.',fontsize=8)
fig.text(.06,.032,'Discovery candidates were not scored by this final benchmark.',fontsize=8)
finish(fig,'Figure_S1_electrostatics_controls')

# Chart contract: paired dot comparison of descriptive RMSDs, no invented passing threshold.
fig,ax=plt.subplots(figsize=(9,5.4))
heading(fig,'S2  TCRmodel2 calibration geometry','Cα RMSD against experimental references · three calibration runs · lower values indicate closer recovery')
ax.set_position([.29,.26,.64,.49])
for i,r in enumerate(tcr):
 a=float(r['pMHC_CA_RMSD_A']);b=float(r['TCR_placement_CA_RMSD_A'])
 ax.plot([a,b],[i,i],color=GRAY,lw=2,zorder=1)
 ax.scatter(a,i,color=BLUE,marker='o',s=55,zorder=3,label='pMHC RMSD' if i==0 else None)
 ax.scatter(b,i,color=GOLD,marker='s',s=50,zorder=3,label='TCR placement RMSD' if i==0 else None)
 ax.annotate(f'{a:.2f}',(a,i),xytext=(0,10),textcoords='offset points',ha='center',fontsize=9,color=BLUE)
 ax.annotate(f'{b:.2f}',(b,i),xytext=(0,-18),textcoords='offset points',ha='center',fontsize=9,color=GOLD)
ax.set_yticks(range(3),['1YMM\nEarlier run','1ZGL\nCorrected exclusion run','2WBJ\nCorrected exclusion run']);ax.invert_yaxis();ax.set_ylim(2.6,-.6)
ax.set_xlim(0,58);ax.set_xlabel('Cα RMSD (Å)');ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
ax.legend(loc='upper left',bbox_to_anchor=(0,1.28),frameon=False,ncol=2,fontsize=9)
fig.text(.06,.11,'Reported coordinate-comparison metrics; no new pass/fail threshold applied. Runs are not biological replicates.',fontsize=8)
fig.text(.06,.065,'Corrected 1ZGL and 2WBJ runs exclude their own pMHC reference templates; 1YMM is the earlier run.',fontsize=8)
finish(fig,'Figure_S2_TCR_docking_controls')

with PdfPages(OUT/'All_panels.pdf') as pdf:
 for name,fig in figures:pdf.savefig(fig)
thumbs=[]
for name,fig in figures:
 im=Image.open(OUT/(name+'.png')).convert('RGB');im.thumbnail((1000,750))
 canvas=Image.new('RGB',(1040,800),'#e5e8eb');canvas.paste(im,((1040-im.width)//2,20));ImageDraw.Draw(canvas).text((20,775),name,fill=INK);thumbs.append(canvas)
sheet=Image.new('RGB',(2080,2400),'white')
for i,im in enumerate(thumbs):sheet.paste(im,((i%2)*1040,(i//2)*800))
sheet.save(OUT/'Panel_overview.png')
contracts={
 'Figure_1':'Workflow diagram; separate screen and targeted review, no proportional funnel; source track_definition.',
 'Figure_2':'Four histograms, 1600 comparisons each, shared bins/axes, source primary_score, focal source-ranked candidates, blue dashed markers.',
 'Figure_3':'Eight-row categorical evidence table; source booleans and arm counts, no composite score, missing evidence distinct from biological negatives.',
 'Figure_S1':'Thirty source gate results in ten-by-three matrix, ranks plus categorical status, missing rendered not evaluable.',
 'Figure_S2':'Three source calibration runs, pMHC and TCR placement RMSD paired dots, shared zero-based axis; no invented threshold.',
 'delivery':'Standalone manuscript panels, vector PDF and editable-text SVG, 400 dpi PNG; white background, blue/gold plus neutrals, labels and shapes support non-color reading.',
}
(OUT/'chart_contracts.json').write_text(json.dumps(contracts,indent=2))
with (OUT/'panel_source_map.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['panel','source_ids','caption_status'])
 for name,ids in [('Figure_1_workflow','track_definition;sequence_rankings;recommendations'),('Figure_2_sequence_rankings','sequence_rankings;recommendations'),('Figure_3_candidate_evidence','candidate_evidence;recommendations'),('Figure_S1_electrostatics_controls','electrostatic_gate;electrostatic_readme'),('Figure_S2_TCR_docking_controls','tcr_metrics;tcr_analysis')]:w.writerow([name,ids,'User to write'])
print('Exported',len(figures),'panels to',OUT)
