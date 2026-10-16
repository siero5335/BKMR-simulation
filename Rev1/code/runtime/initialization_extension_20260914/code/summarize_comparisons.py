# Extracted numerical aggregation only; no draft manuscript generation.
from pathlib import Path
import csv,json,math,statistics,subprocess,re,shutil,zipfile,hashlib,time
import tempfile
out=Path('initialization_extension_20260914')
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def num(x):
    try:return float(x)
    except (TypeError,ValueError):return math.nan
def yes(x):return str(x).lower()=='true'
def median(x):
    x=[v for v in x if math.isfinite(v)];return statistics.median(x) if x else math.nan
def fmt(x,d=3):return f'{x:.{d}f}' if math.isfinite(x) else 'NA'
def json_finite(x):
    if isinstance(x,dict):return {k:json_finite(v) for k,v in x.items()}
    if isinstance(x,list):return [json_finite(v) for v in x]
    if isinstance(x,float) and not math.isfinite(x):return None
    return x
def write(name,rows):
    with (out/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fits=read(out/'fit_summaries.csv');cases=read(out/'cases.csv');status=read(out/'chain_status.csv')
base=[r for r in fits if r['window']=='fixed_5000']
assert len(base)==144 and len(cases)==12 and len(status)==144
sets=('original','alternative_1','alternative_2');budgets=(10000,20000,40000,100000)
lookup={(int(r['task_index']),r['initial_set'],int(r['budget']),r['window']):r for r in fits}
init=[];summary=[]
for c in cases:
    tid,n=int(c['task_index']),int(c['n'])
    for b in budgets:
        g=[lookup[tid,s,b,'fixed_5000'] for s in sets];e=[r for r in g if yes(r['evaluable'])]
        complete=len(e)==3
        means=[num(r['beta_mean']) for r in e]
        init.append(dict(task_index=tid,n=n,budget=b,evaluable_sets=len(e),all_three_evaluable=complete,
            primary_pass_sets=sum(yes(r['primary_pass']) for r in e),beta_pass_sets=sum(yes(r['beta_pass']) for r in e),
            beta_mean_range=max(means)-min(means) if complete else math.nan,
            beta_mean_min=min(means) if complete else math.nan,beta_mean_max=max(means) if complete else math.nan,
            maximum_absolute_beta_error=max(num(r['beta_abs_error']) for r in e) if complete else math.nan,
            differs_in_zero_inclusion=len({yes(r['beta_contains_zero']) for r in e})>1 if complete else None,
            differs_in_truth_inclusion=len({yes(r['beta_contains_truth']) for r in e})>1 if complete else None))
for n in (40,80,160,320):
    for b in budgets:
        g=[r for r in base if int(r['n'])==n and int(r['budget'])==b];e=[r for r in g if yes(r['evaluable'])]
        summary.append(dict(n=n,budget=b,configurations=len(g),evaluable=len(e),unavailable=len(g)-len(e),
            primary_pass=sum(yes(r['primary_pass']) for r in e),beta_pass=sum(yes(r['beta_pass']) for r in e),
            beta_mcse_sd_median=median([num(r['beta_mcse_sd_ratio']) for r in e])))
pair=[]
for c in cases:
    tid=int(c['task_index'])
    for s in sets:
        x,y=lookup[tid,s,10000,'fixed_5000'],lookup[tid,s,100000,'fixed_5000']
        if yes(x['evaluable']) and yes(y['evaluable']):
            pair.append(dict(task_index=tid,n=int(c['n']),initial_set=s,primary_10k=yes(x['primary_pass']),primary_100k=yes(y['primary_pass']),
             beta_mean_10k=num(x['beta_mean']),beta_mean_100k=num(y['beta_mean']),beta_mean_change=num(y['beta_mean'])-num(x['beta_mean']),
             beta_abs_error_10k=num(x['beta_abs_error']),beta_abs_error_100k=num(y['beta_abs_error']),
             beta_mcse_sd_10k=num(x['beta_mcse_sd_ratio']),beta_mcse_sd_100k=num(y['beta_mcse_sd_ratio'])))
common_ids={r['task_index'] for r in init if r['budget']==100000 and r['all_three_evaluable']}
ilook={(r['task_index'],r['budget']):r for r in init}
before=[ilook[i,10000]['beta_mean_range'] for i in sorted(common_ids)]
after=[ilook[i,100000]['beta_mean_range'] for i in sorted(common_ids)]
improved=sum(y<x for x,y in zip(before,after));increased=sum(y>x for x,y in zip(before,after))
trans=dict(paired_configurations=len(pair),nonpass_to_pass=sum(not r['primary_10k'] and r['primary_100k'] for r in pair),
 pass_to_nonpass=sum(r['primary_10k'] and not r['primary_100k'] for r in pair),both_pass=sum(r['primary_10k'] and r['primary_100k'] for r in pair),both_nonpass=sum(not r['primary_10k'] and not r['primary_100k'] for r in pair))
warm=[]
for c in cases:
    tid=int(c['task_index'])
    for s in sets:
        x,y=lookup[tid,s,100000,'fixed_5000'],lookup[tid,s,100000,'half_chain']
        if yes(x['evaluable']) and yes(y['evaluable']):warm.append(dict(task_index=tid,n=int(c['n']),initial_set=s,
          primary_fixed=yes(x['primary_pass']),primary_half=yes(y['primary_pass']),beta_mean_fixed=num(x['beta_mean']),beta_mean_half=num(y['beta_mean']),
          beta_mean_difference_half_minus_fixed=num(y['beta_mean'])-num(x['beta_mean'])))
write('initialization_sensitivity.csv',init);write('budget_summary.csv',summary)
if pair:write('paired_extension.csv',pair)
if warm:write('warmup_sensitivity.csv',warm)
errors=[r for r in status if r['status']!='completed']
full=[r for r in base if int(r['budget'])==100000];full_ok=[r for r in full if yes(r['evaluable'])]
pass100=sum(yes(r['primary_pass']) for r in full_ok)
start10=[r for r in base if int(r['budget'])==10000];ok10=[r for r in start10 if yes(r['evaluable'])]
pass10=sum(yes(r['primary_pass']) for r in ok10)
notes=dict(chain_jobs=144,completed_100k_chains=144-len(errors),sampler_error_chains=len(errors),
 configurations=36,evaluable_10k=len(ok10),primary_pass_10k=pass10,evaluable_100k=len(full_ok),primary_pass_100k=pass100,
 paired_transitions=trans,common_three_initialization_datasets=len(common_ids),common_case_ids=sorted(common_ids),
 median_beta_mean_range_10k=fmt(median(before),6),median_beta_mean_range_100k=fmt(median(after),6),
 range_decreased_cases=improved,range_increased_cases=increased,
 warmup_fixed_only=sum(r['primary_fixed'] and not r['primary_half'] for r in warm),
 warmup_half_only=sum(r['primary_half'] and not r['primary_fixed'] for r in warm),
 warmup_both_pass=sum(r['primary_half'] and r['primary_fixed'] for r in warm),
 warmup_neither=sum(not r['primary_half'] and not r['primary_fixed'] for r in warm),
 errors=errors,summary=summary)
(out/'analysis.json').write_text(json.dumps(json_finite(notes),indent=2,allow_nan=False)+'\n')
