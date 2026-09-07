import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path
OUT=Path(__file__).parent/'figures'; OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10.5,'axes.labelsize':9.5,'legend.fontsize':8.5,'figure.dpi':180,'pdf.fonttype':42,'ps.fonttype':42})
def save(fig,name):
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight',pad_inches=.08)
    fig.savefig(OUT/f'{name}.png',bbox_inches='tight',pad_inches=.08,dpi=220)
    plt.close(fig)

# Main Fig 1: five-stage evidence-to-action pathway
fig,ax=plt.subplots(figsize=(7.2,6.1)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
stages=[
 ('1  Human perturbation','84 paired donors; CTHRC1 rose in 82/84','Independent functional and protein-level follow-up','No clinical or in-vivo universality claim'),
 ('2  Saudi-center context','Whole-LV DCM is a different study-specific object','Bridge with cell-resolved/spatial and protein measurements','Not direct replication or Saudi-population validation'),
 ('3  Candidate-state assessment','ASPN: opposite-transition candidate; WNT5A: pathway candidate','Matched independent + joint perturbations and kinetics','Neither candidate is identified as Schnakenberg v'),
 ('4  Spatial and timing evidence','Transcript organization and sampled response windows are observable','Dense kinetics plus calibrated source/ligand measurements','No physical diffusivity or myocardial hard delay'),
 ('5  Dynamical classification','A stability class requires identified ' + r'$A, B, D, \tau$' + ' and geometry','Propagate parameter uncertainty before classification','No cardiac Turing certificate from transcript patterns')]
ys=[.805,.64,.475,.31,.145]
for (title,finding,nextm,boundary),y in zip(stages,ys):
    x=.04; w=.92; h=.135
    box=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.007,rounding_size=.012',linewidth=1,edgecolor='0.3',facecolor='0.985')
    ax.add_patch(box)
    ax.text(x+.018,y+h-.025,title,ha='left',va='center',weight='bold',fontsize=9.1)
    ax.text(x+.020,y+h-.060,'Finding:',weight='bold',fontsize=7.7,va='top')
    ax.text(x+.155,y+h-.060,finding,fontsize=7.5,va='top')
    ax.text(x+.020,y+h-.088,'Next:',weight='bold',fontsize=7.7,va='top')
    ax.text(x+.155,y+h-.088,nextm,fontsize=7.5,va='top')
    ax.text(x+.020,y+h-.116,'Boundary:',weight='bold',fontsize=7.7,va='top')
    ax.text(x+.155,y+h-.116,boundary,fontsize=7.5,va='top')
ax.text(.5,.985,'Evidence is useful when it changes the next measurement and states the claim boundary.',ha='center',va='top',fontsize=9.1,style='italic')
save(fig,'fig1_evidence_to_action')

# Main Fig 2
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.2,3.35)); fig.subplots_adjust(bottom=.25,wspace=.42)
m,lo,hi=.882,.801,.962
ax1.errorbar([m],[0],xerr=[[m-lo],[hi-m]],fmt='o',capsize=4,elinewidth=1.5,markersize=6)
ax1.axvline(0,linewidth=.8,color='.45'); ax1.set_yticks([]); ax1.set_xlim(.72,1.04); ax1.grid(axis='x',alpha=.22)
ax1.set_xlabel('Paired effect (log2 units)'); ax1.set_title('A  CTHRC1 primary endpoint',loc='left',fontsize=9.5)
ax1.text(m,.13,'82/84 donors increased\n95% bootstrap 0.801-0.962',ha='center',va='bottom',fontsize=8.2)
m,lo,hi=.637,.588,.684
ax2.errorbar([m],[0],xerr=[[m-lo],[hi-m]],fmt='o',capsize=4,elinewidth=1.5,markersize=6)
ax2.axvline(0,linewidth=.8,color='.45'); ax2.set_yticks([]); ax2.set_xlim(.50,.73); ax2.grid(axis='x',alpha=.22)
ax2.set_xlabel('Module difference (equal-z units)'); ax2.set_title('B  Seven-gene repair module',loc='left',fontsize=9.5)
ax2.text(m,.13,'81/84 positive\n95% bootstrap 0.588-0.684',ha='center',va='bottom',fontsize=8.2)
for a in (ax1,ax2):
    a.set_ylim(-.16,.34); a.spines['top'].set_visible(False); a.spines['right'].set_visible(False); a.spines['left'].set_visible(False)
fig.text(.5,.06,'Panels use different units by design; their numerical magnitudes are not directly comparable.',ha='center',fontsize=7.9,style='italic')
save(fig,'fig2_human_pair_result')

# Main Fig 3
fig,ax=plt.subplots(figsize=(7.2,4.6)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
boxes=[(.035,.20,.405,.66),(.56,.20,.405,.66)]
titles=['Controlled primary human\ncardiac fibroblasts','KFSHRC whole left-ventricular\nend-stage DCM tissue']
items=[
 ['84 donor-linked basal/TGF-β1 pairs','24 h controlled perturbation','CTHRC1: 82/84 up; mean +0.882 log2','Multiple repair genes also increased'],
 ['2016 source: 5 DCM + 5 non-failing\nwhole-LV hearts','Detailed Results/S1: 1211 genes\n(597 down, 614 up)','Published S1 lookup:\nTGFB1 signed FC -2.01;\nCTHRC1 not thresholded','2009 report: COL1A1 3.08-fold decrease\n(5 DCM + 4 controls)']]
for (x,y,w,h),title,lines in zip(boxes,titles,items):
    b=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.014,rounding_size=.018',linewidth=1,edgecolor='.3',facecolor='.98'); ax.add_patch(b)
    ax.text(x+w/2,y+h-.07,title,ha='center',va='center',weight='bold',fontsize=9.3)
    yy=y+h-.20
    for line in lines:
        ax.text(x+.03,yy,u'• '+line,ha='left',va='top',fontsize=8.0,wrap=True)
        yy -= 0.04*(line.count('\n')+1) + 0.015
ax.add_patch(FancyArrowPatch((.445,.50),(.555,.50),arrowstyle='<->',mutation_scale=13,linewidth=1.1,color='.3'))
ax.text(.5,.91,'Context transfer',ha='center',va='center',fontsize=8.6,weight='bold')
ax.text(.5,.125,'Practical consequence: use discordance to choose the bridge measurement, not to force a replication label.',ha='center',fontsize=8.8,style='italic')
ax.text(.5,.06,'The Riyadh source establishes a Saudi-center case, not nationality or ancestry of all participants.',ha='center',fontsize=8.1)
save(fig,'fig3_saudi_context_transfer')

# Main Fig 4
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.2,3.25),gridspec_kw={'width_ratios':[1.05,1]}); fig.subplots_adjust(bottom=.25,wspace=.32)
ax1.set_ylim(0,1); ax1.set_xlim(-.5,25.5); ax1.set_yticks([]); ax1.set_xlabel('Hours after perturbation'); ax1.hlines(.5,0,24,linewidth=1,color='.35')
pts=[(2,'RNA 2 h\nearliest sampled 4/4',.72),(6,'Ribo 6 h\nearliest sampled 4/4',.27),(24,'84-pair\n24 h endpoint',.72)]
for x,lab,yv in pts:
    ax1.plot(x,.5,'o'); ax1.vlines(x,.5,yv,linewidth=.9); ax1.text(x,yv+(.035 if yv>.5 else -.035),lab,ha='center',va='bottom' if yv>.5 else 'top',fontsize=7.6)
ax1.set_title('A  Sampled timing is not onset',loc='left',fontsize=9.5)
for s in ['top','right','left']: ax1.spines[s].set_visible(False)
ax2.set_ylim(0,1); ax2.set_xlim(-.1,2.35); ax2.set_yticks([]); ax2.set_xlabel('Fixed-delay parameter δ (h)'); ax2.hlines(.54,0,2.1,linewidth=5,alpha=.35); ax2.plot(.116,.54,'s',markersize=6)
ax2.text(.116,.68,'optimum 0.116 h',ha='center',fontsize=7.8); ax2.text(1.05,.39,'reported profile 0-2.1 h',ha='center',fontsize=8); ax2.set_title('B  Zero remains in the original diagnostic',loc='left',fontsize=9.5)
for s in ['top','right','left']: ax2.spines[s].set_visible(False)
fig.text(.5,.06,'Since δ≥0 puts zero on a boundary, the original profile is used only as a diagnostic; it is not a boundary-calibrated confidence interval.',ha='center',fontsize=7.5,style='italic')
save(fig,'fig4_timing_boundary')

# Supplementary figures
fig,ax=plt.subplots(figsize=(6.4,3.6)); tau=np.array([1,20,40]); reported=np.array([10.360,9.962,9.842]); reproduced=np.array([10.360024,9.962434,9.852307])
ax.plot(tau,reported,'o-',label='Reported'); ax.plot(tau,reproduced,'s--',label='Reproduced')
for x,yv in zip(tau,reported): ax.annotate(f'{yv:.3f}',(x,yv),xytext=(0,7),textcoords='offset points',ha='center',fontsize=8)
ax.set_xlabel(r'Benchmark delay $\tau_A$ (benchmark units)'); ax.set_ylabel('Hopf threshold $b_c$'); ax.set_title('Formulation-specific Alfifi Appendix/Galerkin benchmark'); ax.legend(frameon=False); ax.grid(alpha=.22); ax.margins(x=.08,y=.16)
save(fig,'figS1_alfifi_benchmark')

fig,axes=plt.subplots(1,3,figsize=(7.2,3.05)); fig.subplots_adjust(bottom=.25,wspace=.35)
axes[0].bar(['GSE96991'],[.910],width=.55); axes[0].set_ylim(0,1.2); axes[0].set_ylabel('Native log2 effect'); axes[0].text(0,.96,'3/3 up',ha='center',fontsize=8); axes[0].set_title('Mouse paired lines',fontsize=9)
axes[1].errorbar([0],[.882],yerr=[[.081],[.080]],fmt='o',capsize=3); axes[1].set_xlim(-.6,.6); axes[1].set_xticks([0],['GSE97358']); axes[1].set_ylim(0,1.2); axes[1].text(0,1.02,'82/84 up',ha='center',fontsize=8); axes[1].set_title('84 human donor pairs',fontsize=9)
vals=[1.512,1.354,1.389]; axes[2].plot([24,48,72],vals,'o-'); axes[2].set_ylim(1.2,1.6); axes[2].set_xticks([24,48,72]); axes[2].set_xlabel('hours'); axes[2].set_title('Defined human culture',fontsize=9)
for a in axes: a.grid(axis='y',alpha=.2); a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
fig.text(.5,.07,'Each panel uses its study-specific scale and unit; effect magnitudes are not pooled.',ha='center',fontsize=7.8,style='italic')
save(fig,'figS2_perturbation_support')

fig,ax=plt.subplots(figsize=(6.3,3.25)); sec=['Control','Post-MI 1','Post-MI 2','Post-MI 3']; vals=[-.02,.24,.17,.26]; qs=[.895,.013,.013,.013]; x=np.arange(4); ax.axhline(0,color='.4',linewidth=.8); ax.scatter(x,vals,s=45)
for i,(v,q) in enumerate(zip(vals,qs)): ax.text(i,v+(.03 if v>=0 else -.035),f'I={v:.2f}\nq={q:.3f}',ha='center',va='bottom' if v>=0 else 'top',fontsize=8)
ax.set_xticks(x,sec); ax.set_ylabel("Moran's I for transcript axis"); ax.set_ylim(-.12,.37); ax.grid(axis='y',alpha=.2); ax.set_title('Spatial transcript organization is observable; physical diffusion is not')
for s in ['top','right']: ax.spines[s].set_visible(False)
save(fig,'figS3_spatial_axis')

fig,(a1,a2)=plt.subplots(1,2,figsize=(7.0,3.2)); fig.subplots_adjust(bottom=.25,wspace=.35)
a1.bar(['MI LV vs healthy'],[.926],width=.55); a1.set_ylim(0,1.2); a1.set_ylabel('Human log2CPM difference'); a1.text(0,1.0,'q=0.0288',ha='center',fontsize=8); a1.set_title('Human contrast',fontsize=9.5)
a2.plot([8,60,180],[5.569,5.275,3.709],'o-'); a2.set_xlabel('days post-infarction'); a2.set_ylabel('Pig IZ-RZ log2CPM difference'); a2.set_title('Pig regional contrast',fontsize=9.5)
for a in (a1,a2): a.grid(axis='y',alpha=.2); a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
fig.text(.5,.07,'Separate axes prevent a visual claim of conserved cross-species magnitude.',ha='center',fontsize=7.8,style='italic')
save(fig,'figS4_disease_direction')
