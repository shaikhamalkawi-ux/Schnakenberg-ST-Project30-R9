#!/usr/bin/env python3
"""Independent fail-closed replay of Project 30 four GEO accession pilot.
No user account, no network, no journal manuscript modification.
Requires original gz files in --sources (defaults to source_data).
Does not treat group differences as clinical treatment efficacy, delay, or mechanistic proof.
"""
import argparse,csv,gzip,hashlib,io,itertools,json,math,sys,zipfile
from collections import Counter
from pathlib import Path
from statistics import mean
import numpy as np
from scipy.stats import spearmanr

GENE_TO_ID={'CTHRC1':'ENSG00000164932','ASPN':'ENSG00000106819', 'WNT5A':'ENSG00000114251',
            'COL1A1':'ENSG00000108821','POSTN':'ENSG00000133110','FGFR1':'ENSG00000077782',
            'FSTL3':'ENSG00000070404','MMP2':'ENSG00000087245'}
GENES=tuple(GENE_TO_ID)
INPUTS={
  'GSE202714':('GSE202714_mmoita_normCounts_DESeq2.txt.gz','GSE202714_mmoita_RawCounts_DESeq2.txt.gz'),
  'GSE245825':('GSE245825_DESeq2_normalized_table.txt.gz','GSE245825_raw_counts.txt.gz'),
  'GSE246298':('GSE246298_DESeq2_normalized_table.txt.gz','GSE246298_raw_counts.txt.gz'),
  'GSE234715':('GSE234715_20200114A_readCount.xls.gz',),
}
EXPECTED_SHA={
 'GSE202714_mmoita_normCounts_DESeq2.txt.gz':'903ec77d7d413f1203043ac2070f779ae3e9a6dac8c6e175845b94f227a2419d',
 'GSE202714_mmoita_RawCounts_DESeq2.txt.gz':'70c0071b35e5d9d1e0fe84079f053b1ad0951bc735ec8c7e9c92a2176b68db01',
 'GSE245825_DESeq2_normalized_table.txt.gz':'94f816ba5d343ba49f5423f8c5465462074b9dd1f91014d29772904dddee1ccb',
 'GSE245825_raw_counts.txt.gz':'f6f6ca5d53534fb32c3b34d8f2179701c4e23191ef7e0b2d5fa1f23e9928d0e8',
 'GSE246298_DESeq2_normalized_table.txt.gz':'a423b87fd6874061936174543128dcef05c9bceb368842f94a600a2eb6c37e4e',
 'GSE246298_raw_counts.txt.gz':'db57ba44f1c455f20add50c4b6c73acf34902eda3becd2353ee0f59053076c18',
 'GSE234715_20200114A_readCount.xls.gz':'052dd10acf79ae83788ac6c5f0309b8357b3a5460072814a488837c304c80e32',
}


def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def read_selected(path,geneids=False,full_lib=False):
 wanted=set(GENE_TO_ID.values() if geneids else GENES)
 found={};count=0; lib=None;invalid=0
 with gzip.open(path,'rt',encoding='utf-8-sig',errors='strict',newline='') as f:
  reader=csv.reader(f,delimiter='\t');header=next(reader)
  if geneids:
   assert header[0].strip('"') in ('ID',), header[0]
   cols=header[1:]
  elif header[0] in ('Geneid',''):cols=header[1:]
  else:cols=header  # header is unlabelled for gene symbol; data rows have one more field
  if len(cols)!=len(set(cols)):raise ValueError('duplicated columns')
  if full_lib:lib=np.zeros(len(cols),dtype=np.float64)
  for row in reader:
   count+=1
   if len(row)!=len(cols)+1:raise ValueError(f'Row {count} length {len(row)} vs cols {len(cols)}')
   key=row[0].split('.')[0].strip('"')
   if full_lib:
    v=np.array(row[1:],dtype=float)
    if not np.isfinite(v).all() or (v<0).any():invalid+=1
    else:lib+=v
   if key in wanted:
    if key in found:raise ValueError(f'Duplicate gene row: {key}')
    val=np.array(row[1:],dtype=float)
    if not np.isfinite(val).all() or (val<0).any():raise ValueError('Bad values for '+key)
    found[key]=val
 if full_lib and invalid:raise ValueError(f'Full library had {invalid} invalid rows')
 if full_lib and (lib<=0).any():raise ValueError('A library size is nonpositive')
 if geneids: found={g:found[x] for g,x in GENE_TO_ID.items() if x in found}
 return cols,found,count,lib

def export(path,rows,headers=None):
 if not rows:return
 if headers is None:headers=list(rows[0])
 with open(path,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=headers);w.writeheader();w.writerows(rows)

def perm2(x,y):
 v=np.r_[x,y];n=len(x); obs=abs(np.mean(x)-np.mean(y));diffs=[]
 for k in itertools.combinations(range(len(v)),n):
  xx=v[list(k)];yy=np.delete(v,list(k));diffs.append(abs(np.mean(xx)-np.mean(yy)))
 return sum(q>=obs-1e-12 for q in diffs)/len(diffs)

def bh(p):
 p=np.asarray(p,dtype=float);m=len(p);order=np.argsort(p)
 vals=np.minimum.accumulate((p[order]*m/np.arange(1,m+1))[::-1])[::-1]
 q=np.empty(m);q[order]=np.minimum(1,vals)
 return q

def check_202_names(cols):
 import re
 records=[]
 for col in cols:
  # Normalize only the corrupted printed beta glyph; sample type, treatment and replicate must remain explicit.
  mm=re.fullmatch(r'(?P<source>hiPSC-CF|aHCF|fHCF),\s*(?P<condition>control|TGF-.{1,3}1),\s*replicate (?P<replicate>[123])',col)
  if not mm:raise ValueError('Nonadmissible GSE202714 sample header: '+col)
  records.append((mm['source'],'control' if mm['condition']=='control' else 'tgf',int(mm['replicate'])))
 expected={(src,cond,rep) for src in ('hiPSC-CF','aHCF','fHCF') for cond in ('control','tgf') for rep in (1,2,3)}
 if set(records)!=expected or len(set(records))!=18:raise ValueError('nonunique/missing GSE202714 design')
 return records

def primary202(src,out):
 path=src/INPUTS['GSE202714'][0]; cols,mat,nrows,_=read_selected(path,geneids=True)
 rawcols,rawmat,_,_=read_selected(src/INPUTS['GSE202714'][1],geneids=True)
 if len(cols)!=18 or set(mat)!=set(GENES) or set(rawmat)!=set(GENES):raise ValueError('GSE202714 schema mismatch')
 recs=check_202_names(cols);raw_recs=check_202_names(rawcols)
 if recs!=raw_recs:raise ValueError('GSE202714 raw vs normalized sample ordering mismatch')
 rows=[];vectors=[]
 for gene in GENES:
  vec=np.log2(mat[gene]+1);rv=rawmat[gene]
  for cell in ('hiPSC-CF','aHCF','fHCF'):
   idxc=[i for i,r in enumerate(recs) if r[:2]==(cell,'control')]
   idxt=[i for i,r in enumerate(recs) if r[:2]==(cell,'tgf')]
   x=vec[idxc];y=vec[idxt]; d=float(y.mean()-x.mean())
   rows.append(dict(gene=gene,cell_source=cell,n_control=3,n_treated=3,
        delta_log2_normalized_plus_1=d,mean_control_log2_normp1=float(x.mean()),
        mean_treated_log2_normp1=float(y.mean()),exact_p_two_sided=perm2(x,y),
        n_treated_greater_than_max_control=int((y>x.max()).sum()),
        interpretation='unpaired culture replicate comparison; not donor-level replication'))
  for col,row,val,raw in zip(cols,recs,vec,rv):vectors.append(dict(gene=gene,source='GSE202714',sample_column=col,
      cell_source=row[0],condition=row[1],replicate=row[2],log2_normalized_plus_1=float(val),raw_count=float(raw)))
 for r,q in zip(rows,bh([r['exact_p_two_sided'] for r in rows])):r['q_BH_24_tests']=float(q)
 export(out/'GSE202714_Culture_Source_Contrasts.csv',rows)
 export(out/'GSE202714_Gene_Values.csv',vectors)
 return rows,dict(status='EXPLORATORY_CULTURE_EFFECTS_COMPUTED',n_raw_gene_rows=nrows,n_expression_columns=len(cols),n_genes=len(mat))

def primary245(src,out):
 p=src/INPUTS['GSE245825'][0];cols,mat,nrows,_=read_selected(p)
 rawcols,rawmat,_,_=read_selected(src/INPUTS['GSE245825'][1])
 expected=[f'DCM_{i}' for i in range(1,59)]+[f'Ctrl_{i}' for i in range(1,13)]
 if cols!=expected or rawcols!=expected:raise ValueError('GSE245825 unexpected 58+12 layout')
 if not set(GENES).issubset(mat) or not set(GENES).issubset(rawmat):raise ValueError('Missing source target')
 vals={k:np.log2(v[:58]+1) for k,v in mat.items()};records=[];detail=[]
 for g in GENES:
  for i,(v,n) in enumerate(zip(mat[g],cols)):
   detail.append(dict(gene=g,sample=n,source_dataset='GSE245825' if i<58 else 'GSE116250 reanalysed',
                     group='DCM' if i<58 else 'historical_nonfailing',value_log2_normp1=float(np.log2(v+1))))
 for g in ('CTHRC1','ASPN','WNT5A'):
  for ref in ('COL1A1','POSTN','FSTL3'):
   if g==ref:continue
   rr,pp=spearmanr(vals[g],vals[ref]);records.append(dict(gene=g,comparator=ref,n_DCM=58,
               spearman_within_DCM=float(rr),p_unadjusted=float(pp),
               interpretation='DCM coexpression only; neither a gene-vs-histology relation nor a cross-site control comparison'))
 for r,q in zip(records,bh([r['p_unadjusted'] for r in records])):r['q_BH_'+str(len(records))+'_tests']=float(q)
 export(out/'GSE245825_DCM_WithinCohort_Coexpression.csv',records)
 export(out/'GSE245825_Gene_Values.csv',detail)
 (out/'GSE245825_Histology_HOLD.txt').write_text(
 'HOLD: No verified individual biopsy fibrosis percentages keyed to DCM_1–58/Hh# in downloaded datasets.\n'
 'All 12 historical controls come from GSE116250 and non-failing left ventricle, vs 58 DCM right ventricular septum.\n'
 'Do not treat between-control differences as biological replication; do not claim histology correlation.\n',encoding='utf-8')
 return records,dict(status='COEXPRESSION_COMPUTED_HISTOLOGY_HOLD',n_rows=nrows,n_DCM=58,n_nonfailing=12)

def primary246(src,out):
 p=src/INPUTS['GSE246298'][0];cols,mat,nrows,_=read_selected(p)
 expected=[f'SF{n}-{j}' for n in ('520','521','558') for j in (1,2,3)]
 if cols!=expected or not set(GENES).issubset(mat):raise ValueError('GSE246298 schema')
 recs=[];summary=[]
 for gene in GENES:
  v=np.log2(mat[gene]+1)
  for line in ('SF520','SF521','SF558'):
   a,b,c=[v[cols.index(f'{line}-{j}')] for j in (1,2,3)]
   recs.extend([dict(gene=gene,source_line=line,comparison='fibrotic_vehicle_minus_control',delta_log2normp1=float(b-a)),
                dict(gene=gene,source_line=line,comparison='inhibitor_minus_fibrotic_vehicle',delta_log2normp1=float(c-b))])
 for gene in GENES:
  for comp in ('fibrotic_vehicle_minus_control','inhibitor_minus_fibrotic_vehicle'):
   xx=[r['delta_log2normp1'] for r in recs if r['gene']==gene and r['comparison']==comp]
   summary.append(dict(gene=gene,comparison=comp,n_source_lines=len(xx),mean_delta=float(np.mean(xx)),
           min_delta=float(np.min(xx)),max_delta=float(np.max(xx)),n_positive=sum(x>0 for x in xx),
           interpretation='descriptive 3-source-line contrast; no clinical efficacy or CTHRC1-specific perturbation'))
 export(out/'GSE246298_Organoid_BySource.csv',recs)
 export(out/'GSE246298_Organoid_Summary.csv',summary)
 return summary,dict(status='EXPLORATORY_ORGANOID_EFFECTS_REPRODUCED',n_rows=nrows,n_samples=len(cols))

def primary234(src,out):
 p=src/INPUTS['GSE234715'][0];cols,mat,nrows,full_lib=read_selected(p,full_lib=True)
 expect=[f'P3{x}F{day}s{j}' for day in (2,4) for x in ('','d') for j in (1,2,3)]
 if cols!=expect:raise ValueError('GSE234715 unexpected sample headers')
 # Strict: donor/GSM sample joins from GEO do not prove what "d" in P3dF means.
 # We report gene abundance in each native header without taking case/control contrasts.
 values=[]
 for g in GENES:
  for i,c in enumerate(cols):
   count=mat.get(g,np.zeros(len(cols)))[i]
   v=np.log2((count / full_lib[i]) * 1e6 + 1)
   values.append(dict(gene=g,native_sample_column=c,day_label_from_column=int(c.split('F')[1].split('s')[0]),
      source_group_code=c.split('F')[0].replace('P3',''),log2CPM_plus_1=float(v),raw_count=int(count),
      full_gene_library_size=int(full_lib[i]),GSM_mapping_status='HOLD_NO_OFFICIAL_NATIVE_COL_TO_GSM_JOIN'))
 export(out/'GSE234715_GeneValues_NO_GROUP_COMPARISON.csv',values)
 (out/'GSE234715_Mapping_HOLD.txt').write_text(
 'HOLD: Matrix groups P3F2/P3dF2/P3F4/P3dF4 are not publicly joined to the 12 GEO GSM identifiers in the admitted evidence.\n'
 'Native column data and full gene library normalized log2(CPM+1) extracted, but no assumed F vs dF biological labels.\n',encoding='utf-8')
 return values,dict(status='EXPRESSION_NORMALIZED_ONLY__INJURY_LABELS_HOLD',n_samples=12,n_rows=nrows)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--sources',type=Path,default=Path('source_data'))
 ap.add_argument('--out',type=Path,default=Path('/mnt/data/Project30_GEO_R2_Independent_20261017'))
 args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 verified=[]
 for dataset,fnames in INPUTS.items():
  for fname in fnames:
   path=args.sources/fname
   if not path.exists():raise SystemExit('SOURCE MISSING '+str(path))
   sh=sha(path)
   if sh!=EXPECTED_SHA[fname]:raise SystemExit('SOURCE SHA MISMATCH '+fname)
   with gzip.open(path,'rb') as h:
    while h.read(1<<20):pass
   verified.append(dict(dataset=dataset,file=fname,bytes=path.stat().st_size,sha256=sh,source=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{dataset[:-3]}nnn/{dataset}/suppl/{fname}'))
 export(args.out/'Source_SHA256_Audit.csv',verified)
 with open(args.out/'Ensembl_NCBI_Target_Map.csv','w') as f:
  w=csv.writer(f);w.writerow(['symbol','stable_ensembl_id','reference'])
  for gene,eid in GENE_TO_ID.items():w.writerow([gene,eid,f'https://www.ncbi.nlm.nih.gov/gene/?term={eid}'])
 reports=[]
 for dataset,fun in (('GSE246298',primary246),('GSE202714',primary202),('GSE245825',primary245),('GSE234715',primary234)):
  rr,diag=fun(args.sources,args.out);diag['dataset']=dataset;reports.append(diag)
 (args.out/'RESULTS_REGISTER.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2))
 # Stable check: known original Colab 246298 CTHRC1 positive vehicle and negative inhibitor per 3 sources.
 with open(args.out/'GSE246298_Organoid_Summary.csv') as f:
  rows=list(csv.DictReader(f));one=[x for x in rows if x['gene']=='CTHRC1']
  assert len(one)==2 and abs(float(one[0]['mean_delta'])-1.1886978289578805)<1e-12
  assert abs(float(one[1]['mean_delta'])+1.2330769881611643)<1e-12
 with open(args.out/'GSE202714_Culture_Source_Contrasts.csv') as f:
  rows=list(csv.DictReader(f));assert len(rows)==24 and len([x for x in rows if x['gene']=='CTHRC1'])==3
 with open(args.out/'GSE245825_DCM_WithinCohort_Coexpression.csv') as f:
  rows=list(csv.DictReader(f));assert len(rows)==9
  assert all('q_BH_9_tests' in r for r in rows)
  target=next(r for r in rows if r['gene']=='CTHRC1' and r['comparator']=='COL1A1')
  assert abs(float(target['q_BH_9_tests'])-3.455467599044716e-05)<1e-13
 with open(args.out/'GSE234715_GeneValues_NO_GROUP_COMPARISON.csv') as f:
  rows=list(csv.DictReader(f));assert len(rows)==96
 print('R2 REPLAY COMPLETED',json.dumps(reports,indent=2))

if __name__=='__main__':main()
