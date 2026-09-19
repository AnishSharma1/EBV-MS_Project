from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path('/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/manuscript_figures_v2')
COL = {'ink':'#233142','blue':'#2F6B9A','light_blue':'#E6F0F7','gold':'#B87918','light_gold':'#F7EEDB','gray':'#637381','light_gray':'#F2F4F5','red':'#AE4A45','green':'#2D7A62'}

def save(fig, stem):
    for ext in ('png','pdf','svg'):
        fig.savefig(OUT/f'{stem}.{ext}', dpi=400, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def box(ax, x,y,w,h,text, face, edge=COL['ink'], fs=9, weight='normal'):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.015,rounding_size=0.025',facecolor=face,edgecolor=edge,linewidth=1.2)
    ax.add_patch(p)
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,color=COL['ink'],weight=weight,wrap=True)
    return (x+w,y+h/2)

def arrow(ax, a,b, color=COL['gray']):
    ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='-|>',lw=1.4,color=color,shrinkA=3,shrinkB=3))

# Branched workflow
fig,ax=plt.subplots(figsize=(10.2,4.4))
ax.set_xlim(0,10.8); ax.set_ylim(0,6.5); ax.axis('off')
start=box(ax,.35,2.63,1.45,1.05,'6,400\nallele-specific\ncomparisons',COL['light_blue'],fs=10,weight='bold')
ax.text(3.4,6.12,'Initial evidence-review track',fontsize=10,weight='bold',color=COL['blue'],ha='center')
ax.text(6.6,6.12,'Later expanded-cohort track',fontsize=10,weight='bold',color=COL['gold'],ha='center')
upper1=box(ax,2.5,4.52,1.55,.92,'8 selected for\nevidence audit',COL['light_blue'],edge=COL['blue'])
upper2=box(ax,4.8,4.52,1.7,.92,'2 medium priority\n6 holds\n0 high priority',COL['light_blue'],edge=COL['blue'])
lower1=box(ax,2.5,.82,1.55,.92,'49-pair\nexpanded cohort',COL['light_gold'],edge=COL['gold'])
lower2=box(ax,4.8,.82,1.7,.92,'44 with register\npredictor agreement',COL['light_gold'],edge=COL['gold'])
lower3=box(ax,7.2,.82,1.55,.92,'2 with complete\ncomputational criteria',COL['light_gold'],edge=COL['gold'])
lower4=box(ax,9.35,.82,1.1,.92,'7 provisional\ntiers\n42 untiered',COL['light_gold'],edge=COL['gold'],fs=8)
arrow(ax,start,(2.5,4.98),COL['blue']); arrow(ax,start,(2.5,1.28),COL['gold'])
arrow(ax,(4.05,4.98),(4.8,4.98),COL['blue']); arrow(ax,(4.05,1.28),(4.8,1.28),COL['gold']); arrow(ax,(6.5,1.28),(7.2,1.28),COL['gold']); arrow(ax,(8.75,1.28),(9.35,1.28),COL['gold'])
ax.text(5.65,3.45,'These are separate analyses drawn from the same initial comparison universe.',ha='center',va='center',fontsize=8.5,color=COL['gray'],style='italic')
ax.text(8.0,.26,'No high-priority designation: frozen public-evidence endpoint unavailable',ha='center',va='center',fontsize=7.7,color=COL['red'])
ax.plot([2.2,10.6],[4.12,4.12],color='#D7E4EE',lw=0.8); ax.plot([2.2,10.6],[2.25,2.25],color='#EBDAB3',lw=0.8)
save(fig,'figure_1_branched_workflow')

# Rank sensitivity
fig,ax=plt.subplots(figsize=(7.2,4.2))
methods=['TCR-facing\nBLOSUM62','Full-core\nBLOSUM62','TCR-facing\nidentity','Grantham\ndistance']
x=range(4)
series={
    'DRB1*13:03: BALF5$_{627-641}$-TALDO1$_{108-122}$':[13,173,2,15],
    'DRB1*15:01: BALF5$_{627-641}$-TALDO1$_{216-230}$':[14,3,4,33]
}
colors=[COL['blue'],COL['gold']]
for (label,ranks),c in zip(series.items(),colors):
    ax.plot(x,ranks,color=c,marker='o',ms=7,lw=2.0,label=label)
    for xi,yi in zip(x,ranks):
        ax.annotate(str(yi),(xi,yi),xytext=(0,-13 if yi<30 else 8),textcoords='offset points',ha='center',fontsize=8.5,color=c,weight='bold')
ax.set_xticks(list(x),methods)
ax.set_ylabel('Within-allele rank (1 = highest)')
ax.set_yscale('log'); ax.invert_yaxis(); ax.set_ylim(1600,1)
ax.set_yticks([1,3,10,30,100,300,1000,1600]); ax.set_yticklabels(['1','3','10','30','100','300','1,000','1,600'])
ax.grid(axis='y',color='#D9E0E5',lw=.7); ax.spines[['top','right']].set_visible(False)
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.26),frameon=False,fontsize=8,ncol=1,handlelength=2.3)
ax.tick_params(labelsize=9); fig.subplots_adjust(bottom=.31,left=.16,right=.97,top=.97)
save(fig,'figure_4_rank_sensitivity')
