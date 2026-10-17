"""Post-hoc exploratory sensitivity: exact source-stratified 3+3 permutation and DCM LOO rho."""
import csv,itertools,json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
P=Path(__file__).resolve().parent

def readcsv(p):
 with open(p,encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def outcsv(path,rows):
 with open(path,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def bh(p):
 arr=np.asarray(p,float);idx=np.argsort(arr);val=np.minimum.accumulate((arr[idx]*len(arr)/np.arange(1,len(arr)+1))[::-1])[::-1]
 q=np.empty(len(arr));q[idx]=np.clip(val,0,1);return q

a=readcsv(P/'GSE202714_Gene_Values.csv');gene_order=sorted({r['gene'] for r in a});res=[]; all_perm=[]
for gene in gene_order:
 this=[r for r in a if r['gene']==gene]
 perm_per=[];observed=[]; max_dir=[]
 for src in ('hiPSC-CF','aHCF','fHCF'):
  six=np.array([float(r['log2_normalized_plus_1']) for r in this if r['cell_source']==src])
  labels=[r['condition'] for r in this if r['cell_source']==src]
  if len(six)!=6 or labels.count('tgf')!=3 or labels.count('control')!=3:raise ValueError('source condition schema')
  y=six[np.array(labels)=='tgf'];x=six[np.array(labels)=='control'];observed.append(float(y.mean()-x.mean()))
  diffs=[]
  for inds in itertools.combinations(range(6),3):
   m=np.zeros(6,bool);m[list(inds)]=True
   diffs.append(float(six[m].mean()-six[~m].mean()))
  perm_per.append(np.asarray(diffs))
 # 20^3=8,000 balanced assignments preserving each source (3 per group).
 v=(perm_per[0][:,None,None]+perm_per[1][None,:,None]+perm_per[2][None,None,:])/3
 obs=float(np.mean(observed));p=float(np.mean(np.abs(v)>=abs(obs)-1e-12))
 res.append(dict(gene=gene,delta_mean_across_source_types=obs,p_exact_block_stratified_two_sided=p,
  n_treatment_positives=sum(x>0 for x in observed),min_source_effect=min(observed),max_source_effect=max(observed),
  sensitivity_min_2source_mean=min(float(np.mean(x)) for x in itertools.combinations(observed,2)),
  sensitivity_max_2source_mean=max(float(np.mean(x)) for x in itertools.combinations(observed,2)),
  test_scope='Exploratory; treatment-label permutations within source; 3+3 cultures per source; no donor-level inference'))
for r,q in zip(res,bh([r['p_exact_block_stratified_two_sided'] for r in res])):r['q_BH_8_genes']=float(q)
outcsv(P/'GSE202714_ExactSourceStratified_Exploratory.csv',res)

b=readcsv(P/'GSE245825_Gene_Values.csv');groups={}
for r in b:
 if r['group']=='DCM':groups.setdefault(r['gene'],{})[r['sample']]=float(r['value_log2_normp1'])
ref=['COL1A1','POSTN','FSTL3'];targets=['CTHRC1','ASPN','WNT5A']
rr=[]
for gene in targets:
 for comp in ref:
  v1=np.asarray([groups[gene][f'DCM_{i}'] for i in range(1,59)])
  v2=np.asarray([groups[comp][f'DCM_{i}'] for i in range(1,59)])
  r0,p0=spearmanr(v1,v2)
  lo=[spearmanr(np.delete(v1,i),np.delete(v2,i)).statistic for i in range(len(v1))]
  rng=np.random.default_rng(20261017)
  boots=[]
  for _ in range(3000):
   inds=rng.integers(0,58,58)
   boots.append(spearmanr(v1[inds],v2[inds]).statistic)
  rr.append(dict(gene=gene,comparator=comp,n=58,rho=float(r0),p_two_sided=float(p0),
    loo_min_rho=float(min(lo)),loo_max_rho=float(max(lo)),
    bootstrap_percentile_95_lower=float(np.nanquantile(boots,.025)),
    bootstrap_percentile_95_upper=float(np.nanquantile(boots,.975)),
    bootstrap_resamples=3000,seed=20261017,
    interpretation='sample-level within DCM coexpression; independent patient IDs/linked fibrosis unknown'))
for r,q in zip(rr,bh([x['p_two_sided'] for x in rr])):r['q_BH_9_tests']=float(q)
outcsv(P/'GSE245825_Coexpression_LOO_Bootstrap.csv',rr)
print('SOURCE STRATIFIED')
for r in res:print(r['gene'],round(r['delta_mean_across_source_types'],3),'p',r['p_exact_block_stratified_two_sided'],'q',r['q_BH_8_genes'])
print('COEXPRESSION')
for r in rr[:3]:print(r['gene'],r['comparator'],'rho',r['rho'],'95CI',r['bootstrap_percentile_95_lower'],r['bootstrap_percentile_95_upper'])
