import csv, statistics as st, json, pathlib, zipfile
root=pathlib.Path(__file__).resolve().parents[1]
import argparse
parser=argparse.ArgumentParser(description='Reaggregate published CSVs; no MCMC.')
parser.add_argument('--output-dir',type=pathlib.Path,default=root/'_work/coauthor_checks')
dest=parser.parse_args().output_dir.resolve()
if dest==root/'coauthor_additions':raise SystemExit('Use a separate output directory.')
dest.mkdir(parents=True,exist_ok=True)
def read(p):
    with (root/p).open() as f: return list(csv.DictReader(f))
def num(r,k): return float(r[k])
fits=[r for r in read('initialization_extension/fit_summaries.csv') if r['window']=='fixed_5000']
chains=read('initialization_extension/beta_by_chain.csv')
summary={}
for flag in ['TRUE','FALSE']:
    rows=[r for r in fits if r['primary_pass']==flag]
    ranges=[]
    for r in rows:
        cc=[num(c,'beta_mean') for c in chains if all(c[k]==r[k] for k in ['task_index','initial_set','budget','window'])]
        v=(max(cc)-min(cc))/num(r,'beta_sd')
        assert abs(v-num(r,'chain_mean_range_over_sd'))<1e-10
        ranges.append(v)
    summary['chain_range_'+flag]={'n':len(rows),'median':st.median(ranges),'max':max(ranges)}
refs=[]
for r in fits:
    if r['budget']=='100000': continue
    rr=[s for s in fits if s['task_index']==r['task_index'] and s['initial_set']!=r['initial_set'] and s['budget']=='100000']
    if len(rr)!=2 or any(s['primary_pass']!='TRUE' for s in rr):continue
    refs.append(dict(task=r['task_index'],initial_set=r['initial_set'],budget=r['budget'],flag=r['primary_pass'],**{k:abs(num(r,k)-st.mean(num(s,k) for s in rr))/num(r,'beta_sd') for k in ['beta_mean','beta_q975']}))
summary['references']={'n':len(refs),'datasets':len(set(r['task'] for r in refs))}
for flag in ['TRUE','FALSE']:
    rs=[r for r in refs if r['flag']==flag]
    summary['reference_'+flag]={'n':len(rs),**{k:{'median':st.median(r[k] for r in rs),'max':max(r[k] for r in rs)} for k in ['beta_mean','beta_q975']}}
b=read('coauthor_additions/beta_estimand_by_task.csv')
for status in sorted(set(r['primary_status'] for r in b)):
    rr=[r for r in b if r['primary_status']==status]
    summary['endpoint_'+status]={'n':len(rr),'above':sum(max(num(r,'mcse_q025'),num(r,'mcse_q975'))>.1*num(r,'interval_width')/2 for r in rr)}
tt=read('coauthor_additions/task_table.csv')
different=[r for r in tt if r['chain_seeds']!=tt[0]['chain_seeds']]
summary['different_seeds']={'count':len(different),'groups':sorted(set(r['bundle'] for r in different))}
aa=read('coauthor_additions/all_task_results.csv')
summary['matched']={}
for bundle in ['balanced_main','default_support']:
    for n in ['40','80']:
        rr=[r for r in aa if r['bundle']==bundle and r['n']==n and r['iter']=='10000' and r['varsel']==('TRUE' if n=='40' else 'FALSE')]
        summary['matched'][bundle+'_'+n]={'count':len(rr),**{k:st.median(num(r,k) for r in rr) for k in ['max_rhat','min_bulk_ess']},'event_rate':st.median(num(r,'label_sum')/int(n) for r in rr)}
out=json.dumps(summary,indent=2,ensure_ascii=False)
(dest/'numeric_checks.json').write_text(out)
print(out)

def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write_csv(dest/'reference_comparisons.csv',refs)
write_csv(dest/'chain_ranges.csv',[{k:r[k] for k in ['task_index','n','initial_set','budget','primary_pass','beta_sd','chain_mean_range','chain_mean_range_over_sd']} for r in fits])
