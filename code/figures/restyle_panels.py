from pathlib import Path
import csv, json, hashlib, shutil, zipfile, sys
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image, ImageDraw

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'outputs/manuscript_final_panels_2026-09-05'
CLEAN='--clean' in sys.argv
OUT=ROOT/('outputs/manuscript_figures_only_2026-09-06' if CLEAN else 'outputs/manuscript_restyled_panels_2026-09-05')
OUT.mkdir(exist_ok=True)
DATA=OUT/'source_data';DATA.mkdir(exist_ok=True)
for p in (OLD/'source_data').iterdir():
 if p.is_file():shutil.copyfile(p,DATA/p.name)
GATE=Path('/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/processed/pmhc_surface_electrostatics_v2_controls_2026-08-30/control_gate.json')
shutil.copyfile(GATE,DATA/'control_gate.json')
def read(name,delimiter=','):
 with (DATA/name).open() as f:return list(csv.DictReader(f,delimiter=delimiter))
rank=read('sequence_rankings.csv');ev=read('candidate_evidence.csv');gate=read('electrostatic_gate.csv');tcr=read('tcr_metrics.tsv','\t')
cg=json.loads((DATA/'control_gate.json').read_text())
assert len(rank)==6400 and Counter(r['allele'] for r in rank)=={f'HLA-DRB1*{a}':1600 for a in ['03:01','08:01','13:03','15:01']}
assert Counter(r['stage1_status'] for r in ev)=={'stage1_hold':6,'stage1_medium_priority':2}
leads={r['allele']:next(x for x in rank if x['pair_id']==r['pair_id']) for r in ev if r['stage1_status']=='stage1_medium_priority'}
assert [int(leads['HLA-DRB1*'+a]['hla_rank']) for a in ['13:03','15:01']]==[13,14]
BLUE='#286FA6'; DARK='#173D5B'; GOLD='#C37A1F'; INK='#162D3D'; MUTED='#586C7B'; PALE='#E6F0F7'; WARM='#FAEBD6'; GRAY='#EDF0F2'; WHITE='#FFFFFF'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':INK,'axes.edgecolor':MUTED,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
figs=[]
def header(fig,num,title,sub):
 if CLEAN:return
 fig.text(.045,.962,num,fontsize=24,weight='bold',color=BLUE,va='top')
 fig.text(.13,.953,title,fontsize=18,weight='bold',va='top')
 fig.text(.13,.897,sub,fontsize=9.5,color=MUTED,va='top')
 fig.add_artist(Line2D([.045,.955],[.855,.855],transform=fig.transFigure,color=BLUE,lw=1.2))
def note(fig,lines):
 if CLEAN:return
 for i,line in enumerate(lines):fig.text(.045,.084-i*.027,line,fontsize=8.2,color=MUTED)
def export(fig,name):
 if CLEAN:
  for ax in fig.axes:
   title=ax.get_title(loc='left')
   if len(title)>4 and title[0] in 'ABCD' and title[1:4]=='   ':ax.set_title(title[4:],loc='left')
 for ext in ['png','svg','pdf']:fig.savefig(OUT/(name+'.'+ext),dpi=400,**({'bbox_inches':'tight','pad_inches':.10} if CLEAN else {}))
 figs.append((name,fig))

# Contract: branching workflow; large counts for emphasis, diagram geometry does not encode counts.
fig=plt.figure(figsize=(10,7.3));header(fig,'01','Study workflow','Parallel sequence and structural analyses, followed by a targeted evidence review')
ax=fig.add_axes([.04,.135,.92,.68]);ax.set(xlim=(0,10),ylim=(0,7));ax.axis('off')
def tile(x,y,w,h,title,body,kind='blue'):
 color=BLUE if kind=='blue' else GOLD
 ax.add_patch(Rectangle((x,y),w,h,facecolor=PALE if kind=='blue' else WARM,edgecolor='none'))
 ax.add_patch(Rectangle((x,y+h-.10),w,.10,facecolor=color,edgecolor='none'))
 ax.text(x+.16,y+h-.26,title,weight='bold',fontsize=11,color=color,va='top')
 ax.text(x+.16,y+.18,body,fontsize=9.5,va='bottom',linespacing=1.4)
def arrow(a,b,dash=False):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,color=MUTED,lw=1.2,linestyle='--' if dash else '-'))
tile(1.9,5.95,6.2,.95,'Peptide library','IEDB records + canonical human protein sequences')
tile(.1,3.95,4.35,1.3,'6,400 sequence comparisons','1,600 per allele\nDRB1*03:01 · *08:01 · *13:03 · *15:01')
tile(5.55,3.95,4.35,1.3,'Structural analysis','Predicted peptide–HLA complexes\nSeparate V3 structural shortlist')
arrow((3.4,5.95),(2.27,5.25));arrow((6.6,5.95),(7.72,5.25))
tile(.1,1.95,4.35,1.25,'8 candidates reviewed','Binding · registers · provenance\nAdditional evidence and gaps')
tile(5.55,1.95,4.35,1.25,'Control benchmarks','Surface electrostatics\nTCR docking')
arrow((2.27,3.95),(2.27,3.2),True);ax.text(2.45,3.57,'Targeted follow-up',fontsize=8.5,va='center',color=MUTED)
arrow((7.72,3.95),(7.72,3.2))
tile(.1,.03,2.15,1.1,'2 medium priority','Binding / register\nstudies proposed')
tile(2.55,.03,1.90,1.1,'6 on hold','Evidence gaps',kind='gold')
tile(5.55,.03,4.35,1.1,'Control limitations','No additional candidate support',kind='gold')
arrow((1.35,1.95),(1.18,1.13));arrow((3.4,1.95),(3.5,1.13));arrow((7.72,1.95),(7.72,1.13))
note(fig,['Dashed connector indicates a subsequent targeted review, not a unified all-feature ranking.','Experimental peptide–HLA binding, register confirmation and T-cell recognition remain pending.'])
export(fig,'Figure_1_workflow')

# Contract: gallery-inspired classic histograms, identical bins and axes, lead marker separated from text.
fig,axs=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
header(fig,'02','Within-allele sequence scores','6,400 comparisons · predicted TCR-facing positions P2, P3, P5, P7 and P8')
fig.subplots_adjust(left=.095,right=.96,top=.745,bottom=.16,hspace=.54,wspace=.19)
vals=[float(r['primary_score']) for r in rank];bins=np.linspace(min(vals),max(vals),31)
for ax,a,letter in zip(axs.flat,['03:01','08:01','13:03','15:01'],'ABCD'):
 allele='HLA-DRB1*'+a;v=[float(r['primary_score']) for r in rank if r['allele']==allele]
 ax.hist(v,bins=bins,color=BLUE,alpha=.83,edgecolor=WHITE,lw=.7,zorder=2)
 ax.stairs(np.histogram(v,bins)[0],bins,color=DARK,lw=.8,zorder=3)
 ax.set_ylim(0,255);ax.set_xlim(bins[0]-.025,bins[-1]+.025)
 ax.grid(axis='y',color='#DDE5EA',lw=.7);ax.set_axisbelow(True);ax.tick_params(length=3)
 ax.set_title(f'{letter}   HLA-DRB1*{a}',loc='left',fontsize=11,weight='bold',pad=15)
 ax.text(1,1.065,'n = 1,600',transform=ax.transAxes,ha='right',fontsize=9,color=MUTED)
 if allele in leads:
  r=leads[allele];x=float(r['primary_score'])
  ax.axvline(x,color=GOLD,ls=(0,(4,3)),lw=1.7,zorder=4)
  ax.scatter([x],[246],s=25,marker='D',color=GOLD,zorder=5)
  ax.annotate(f"BALF5–TALDO1\nRank {r['hla_rank']}/1,600 · score {x:.3f}",xy=(x,.97),xycoords=ax.get_xaxis_transform(),xytext=(-10,0),textcoords='offset points',ha='right',va='top',fontsize=9,color=INK,fontweight='medium')
for ax in axs[:,0]:ax.set_ylabel('Pair comparisons',labelpad=8)
for ax in axs[1,:]:ax.set_xlabel('BLOSUM62-based similarity score',labelpad=8)
note(fig,['Higher score indicates greater similarity, not percentage identity or probability. Ranks are within each allele.','Geometry and binding percentiles were not used in these original sequence ranks.'])
export(fig,'Figure_2_sequence_rankings')

# Contract: annotated categorical heatmap; four booleans, two count columns, one outcome. No new score.
fig=plt.figure(figsize=(12,7.5));header(fig,'03','Candidate evidence matrix','Eight preselected peptide pairs · source-defined consensus and exact-HLA evidence')
ax=fig.add_axes([.04,.20,.92,.56]);ax.set(xlim=(-3.3,8.4),ylim=(8.5,-1.1));ax.axis('off')
fields=['ebv_binding_consensus','self_binding_consensus','ebv_register_consensus_matches_declared','self_register_consensus_matches_declared','iedb_exact_hla_positive_arm_count','immunopeptidome_exact_hla_compatible_arm_count']
xs=[0,1,2.25,3.25,4.6,5.85]
for a,b,label in [(0,1.92,'Binding consensus'),(2.25,4.17,'Register agreement'),(4.6,6.77,'Evidence arms / 2')]:
 ax.plot([a,b],[-.68,-.68],color=BLUE,lw=2);ax.text((a+b)/2,-.82,label,ha='center',va='bottom',fontsize=10,weight='bold',color=BLUE)
for x,label in zip(xs,['EBV','Self','EBV','Self','IEDB','Ligand']):ax.text(x+.46,-.32,label,ha='center',fontsize=9,color=MUTED)
ax.text(-3.25,-.32,'Peptide pair / HLA-DRB1',fontsize=9,color=MUTED)
ax.text(7.55,-.32,'Stage 1',ha='center',fontsize=9,color=MUTED)
for i,r in enumerate(ev):
 medium=r['stage1_status']=='stage1_medium_priority'
 label=r['ebv_protein'].replace('BALF4_gB','BALF4')+'–'+r['self_protein']
 ax.text(-3.2,i+.26,label,fontsize=10,weight='bold' if medium else 'normal',va='center',color=BLUE if medium else INK)
 ax.text(-3.2,i+.62,r['allele'].replace('HLA-DRB1',''),fontsize=8.5,color=MUTED,va='center')
 for j,(field,x) in enumerate(zip(fields,xs)):
  raw=r[field];positive=raw=='True' if j<4 else int(raw)>0
  ax.add_patch(Rectangle((x,i+.05),.92,.82,facecolor=BLUE if positive else GRAY,edgecolor='white',lw=1))
  ax.text(x+.46,i+.46,('Yes' if positive else 'No') if j<4 else raw,ha='center',va='center',fontsize=10,weight='bold' if positive else 'normal',color='white' if positive else MUTED)
 ax.add_patch(Rectangle((7.02,i+.05),1.2,.82,facecolor=PALE if medium else WARM,edgecolor='none'))
 ax.text(7.62,i+.46,'Medium\npriority' if medium else 'Hold',ha='center',va='center',fontsize=9,weight='bold',color=BLUE if medium else GOLD)
note(fig,['Register agreement: consensus matches the declared core. “No” is lack of consensus, not proven nonbinding.','Evidence counts refer to the two peptide arms. Zero compatible ligand records indicates missing evidence.','Medium priority applies to initial binding and register studies; it does not establish T-cell cross-reactivity.'])
export(fig,'Figure_3_candidate_evidence')

# Contract: annotated status tiles based directly on saved control_gate failures, not primary rank alone.
failure_keys={(r['layer'],r['pair_id'],r['endpoint']) for r in cg['failures']}
lookup={(r['layer'],r['pair_id'],r['endpoint']):r for r in gate}
keys=sorted({(r['layer'],r['pair_id']) for r in gate},key=lambda k:({'pdb':0,'af_271828':1,'af_314159':2}[k[0]],k[1]))
names={'PAIR_HY2E11_BALF5_MBP':'Hy.2E11 · BALF5 / MBP','PAIR_OB1A12_ENGA_MBP':'Ob.1A12 · ENGA / MBP','PAIR_HY1B11_UL15_MBP':'Hy.1B11 · UL15 / MBP','PAIR_HY1B11_PMM_MBP':'Hy.1B11 · PMM / MBP'}
fig=plt.figure(figsize=(10,8.4));header(fig,'S1','Electrostatics control benchmark','Final control-first V2 analysis · overall gate: FAIL · candidates not evaluated')
ax=fig.add_axes([.04,.19,.92,.61]);ax.axis('off');ax.set(xlim=(-3.45,6.15),ylim=(10.5,-.8))
for j,s in enumerate(['Peptide\nelectrostatics','Composite\nelectrostatics','Surface\nshape']):ax.text(j*2+.87,-.25,s,ha='center',va='bottom',fontsize=10,weight='bold',color=BLUE)
actual_status=[]
for i,(layer,pair) in enumerate(keys):
 ax.text(-3.38,i+.28,names[pair],fontsize=9.5,va='center')
 ax.text(-3.38,i+.65,{'pdb':'Experimental PDB','af_271828':'AF seed 271828','af_314159':'AF seed 314159'}[layer],fontsize=8,color=MUTED,va='center')
 for j,end in enumerate(['peptide','composite','shape']):
  k=(layer,pair,end);r=lookup[k];status='NE' if r['status']=='missing' else ('FAIL' if k in failure_keys else 'PASS');actual_status.append(status)
  fill={'NE':GRAY,'FAIL':WARM,'PASS':BLUE}[status];color={'NE':MUTED,'FAIL':GOLD,'PASS':'white'}[status]
  ax.add_patch(Rectangle((j*2,i+.06),1.75,.82,facecolor=fill,edgecolor='white',lw=1.3))
  text='Not evaluable' if status=='NE' else f"Rank {int(r['rank'])}  ·  {status}"
  ax.text(j*2+.875,i+.47,text,ha='center',va='center',color=color,fontsize=9,weight='bold' if status!='NE' else 'normal')
assert Counter(actual_status)=={'FAIL':15,'NE':12,'PASS':3}
note(fig,['Pass/fail includes rank, sensitivity and resampling requirements. Panel sizes differ across controls.','Not evaluable: required comparison-panel coverage was not met. Missing panels are not biological negatives.','Low primary rank alone does not establish a passing control result.'])
export(fig,'Figure_S1_electrostatics_controls')

# Contract: gallery-inspired lollipops, faceted by metric; same zero baseline and shared scale.
fig,axs=plt.subplots(1,2,figsize=(10,6.6),sharex=True,sharey=True)
header(fig,'S2','TCRmodel2 calibration geometry','Three calibration runs · Cα RMSD against experimental references · lower is closer')
fig.subplots_adjust(left=.23,right=.95,top=.72,bottom=.25,wspace=.22)
labels=['1YMM\nEarlier run','1ZGL\nCorrected exclusion','2WBJ\nCorrected exclusion']
for ax,key,title,color,marker in zip(axs,['pMHC_CA_RMSD_A','TCR_placement_CA_RMSD_A'],['A   pMHC geometry','B   TCR placement'],[BLUE,GOLD],['o','D']):
 for i,r in enumerate(tcr):
  val=float(r[key]);ax.hlines(i,0,val,color=color,alpha=.5,lw=3);ax.scatter(val,i,s=90,color=color,marker=marker,zorder=3,edgecolor='white',lw=.8)
  ax.annotate(f'{val:.2f}',(val,i),xytext=(0,13),textcoords='offset points',ha='center',fontsize=11,weight='bold',color=color)
 ax.set_title(title,loc='left',fontsize=11,weight='bold',color=color,pad=20)
 ax.set_xlim(0,58);ax.set_ylim(2.5,-.55);ax.set_yticks(range(3),labels);ax.set_xticks([0,10,20,30,40,50]);ax.set_xlabel('Cα RMSD (Å)',labelpad=10)
 ax.grid(axis='x',color='#DFE6EB',lw=.7);ax.set_axisbelow(True);ax.tick_params(axis='y',length=0)
note(fig,['Descriptive coordinate-comparison metrics; no new pass/fail cutoff. Runs are not biological replicates.','Corrected 1ZGL and 2WBJ runs excluded their own pMHC reference templates; 1YMM is the earlier run.'])
export(fig,'Figure_S2_TCR_docking_controls')

with PdfPages(OUT/'All_panels.pdf') as pdf:
 for name,fig in figs:pdf.savefig(fig,**({'bbox_inches':'tight','pad_inches':.10} if CLEAN else {}))
sheet=Image.new('RGB',(2200,2490),'#E6ECF0')
for i,(name,fig) in enumerate(figs):
 im=Image.open(OUT/(name+'.png')).convert('RGB');im.thumbnail((1060,770))
 x=(i%2)*1100+(1100-im.width)//2;y=(i//2)*830+15
 sheet.paste(im,(x,y))
 if not CLEAN:ImageDraw.Draw(sheet).text(((i%2)*1100+25,(i//2)*830+800),name,fill=INK)
sheet.save(OUT/'Panel_overview.png')
with (OLD/'panel_source_map.csv').open() as f: source_map=list(csv.DictReader(f))
for row in source_map:
 if row['panel']=='Figure_S1_electrostatics_controls':row['source_ids']+=';control_gate'
with (OUT/'panel_source_map.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=list(source_map[0]));writer.writeheader();writer.writerows(source_map)
manifest=[]
for p in sorted(DATA.iterdir()):manifest.append({'snapshot':str(p.relative_to(OUT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
(OUT/'chart_contracts.json').write_text(json.dumps({'reference':'https://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html','renderer':'matplotlib standalone manuscript exports','style':'classic histogram facets, lollipop facets, annotated categorical tiles; blue/gold palette; explicit text and marker shapes','unchanged_data':'identical input score rows, bins, candidate outcomes and RMSD values','correction':'S1 gate statuses now read from source control_gate.json failures; 15 fail, 12 not evaluable, 3 pass','output':'five panels; PDF, editable-text SVG and 400-dpi PNG'},indent=2))
from pypdf import PdfReader
assert len(PdfReader(OUT/'All_panels.pdf').pages)==5
for name,fig in figs:
 assert len(PdfReader(OUT/(name+'.pdf')).pages)==1
 im=Image.open(OUT/(name+'.png'));assert min(im.size)>=(1200 if CLEAN else 2600)
 assert (OUT/(name+'.svg')).stat().st_size>1000
(OUT/'QA.json').write_text(json.dumps({'data_counts':'pass','lead_ranks':'13/1600 and 14/1600','control_gate_status_counts':dict(Counter(actual_status)),'pdf_page_count':5,'panel_formats':['PDF','SVG','PNG'],'visual_review':'Contact sheet inspected; Figure 2 footer simplified to avoid axis-label overlap. Final Figure 2 rendered for inspection.'},indent=2))
with zipfile.ZipFile(OUT/'Restyled_panels.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(OUT.iterdir()):
  if p.suffix in ['.pdf','.svg','.png']:z.write(p,p.name)
print(OUT)
