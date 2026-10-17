#!/usr/bin/env python3
"""Post-hoc, culture-unit sensitivity for already-admitted GSE152250 CTHRC1 values.

Independent timepoint groups: three control and three TGFb1 culture replicates.
For each timepoint, remove one culture from *each* treatment arm (9 combinations);
report range of log2 ratios of arithmetic mean untransformed FPKM.
For the exploratory treatment-by-time contrast (6h to 24h), use a transformed
mean difference per timepoint log2(FPKM+0.5) and exact within-timepoint label
permutations (choose(6,3)^2=400 arrangements). This is a conditional
exchangeability reference, not a claim about biological randomization.
"""
from pathlib import Path
import csv, itertools, math, json, hashlib
from collections import defaultdict
src=Path(__file__).resolve().parent/'results'/'GSE152250_CTHRC1_sample_values.csv'
rows=list(csv.DictReader(src.open()))
assert len(rows)==18
vals=defaultdict(lambda: defaultdict(list))
for r in rows:
    t=int(r['time_h']); arm=r['group']
    assert arm in {'CTRL','TGFb1'} and t in {6,24,48}
    vals[t][arm].append(float(r['CTHRC1_FPKM']))
assert all(len(vals[t][arm])==3 for t in [6,24,48] for arm in ['CTRL','TGFb1'])

def avg(x):return sum(x)/len(x)
def fold_log2(t, ctrl, treated):return math.log2(avg(treated)/avg(ctrl))
def mean_difference(t, ctrl, treated):return avg([math.log2(v+0.5) for v in treated])-avg([math.log2(v+0.5) for v in ctrl])

records=[];perms={}
for t in [6,24,48]:
    ctrl=vals[t]['CTRL'];tr=vals[t]['TGFb1']
    effects=[]
    for i in range(3):
        for j in range(3):
            effects.append(fold_log2(t,[v for k,v in enumerate(ctrl) if k!=i],[v for k,v in enumerate(tr) if k!=j]))
    pooled=ctrl+tr
    perms[t]=[]
    for idxs in itertools.combinations(range(6),3):
        s=set(idxs)
        new_tr=[pooled[i] for i in range(6) if i in s]
        new_ctrl=[pooled[i] for i in range(6) if i not in s]
        perms[t].append(mean_difference(t,new_ctrl,new_tr))
    records.append({'time_h':t,'n_control':3,'n_treated':3,'reference_log2_ratio_of_means':fold_log2(t,ctrl,tr),'leave_one_each_arm_n':len(effects),'loo_min_log2_ratio':min(effects),'loo_max_log2_ratio':max(effects),'loo_all_positive':all(x>0 for x in effects),'mean_difference_of_log2_FPKM_plus_0_5':mean_difference(t,ctrl,tr)})
assert len(perms[6])==len(perms[24])==20
ob=records[1]['mean_difference_of_log2_FPKM_plus_0_5']-records[0]['mean_difference_of_log2_FPKM_plus_0_5']
null=[b-a for a in perms[6] for b in perms[24]]
num=sum(abs(x)>=abs(ob)-1e-12 for x in null)
assert len(null)==400
p=num/len(null)
assert abs(p-0.065)<1e-12
assert records[1]['loo_all_positive'] and records[2]['loo_all_positive'] and not records[0]['loo_all_positive']
result={'classification':'exploratory_posthoc_sensitivity_only','dataset':'GSE152250','input_csv_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'sample_level_units':'18 independent culture records; three per arm at each of 6, 24, 48 h','no_donor_pairing_assumed':True,'no_cross_study_pooling':True,'interaction_test':{'time_comparison':'24h minus 6h','estimand':'difference of arm-mean log2(FPKM+0.5) treatment contrasts','delta':ob,'reference':'conditional within-timepoint 3/3 label exchangeability','enumerated_permutations':400,'absolute_extreme_count':num,'two_sided_exact_P':p,'not_prespecified':True,'not_mechanistic_delay_or_onset':True},'timepoint_sensitivity':records}
out=Path(__file__).resolve().parent/'results'/'GSE152250_Exploratory_Temporal_Sensitivity.json';out.write_text(json.dumps(result,indent=2)+'\n')
rows_out=Path(__file__).resolve().parent/'results'/'GSE152250_Exploratory_LOO.csv'
with rows_out.open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(records[0].keys()));writer.writeheader();writer.writerows(records)
print(json.dumps({'n_t':len(records),'n_arrangements':400,'p':p,'looc_range24':[records[1]['loo_min_log2_ratio'],records[1]['loo_max_log2_ratio']],'looc_range48':[records[2]['loo_min_log2_ratio'],records[2]['loo_max_log2_ratio']]},indent=2))
