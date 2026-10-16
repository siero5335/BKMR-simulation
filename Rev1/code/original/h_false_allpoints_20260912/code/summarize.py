from pathlib import Path
import csv,json,math,statistics

out=Path('h_false_allpoints_20260912')
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def write(name,rows):
    with (out/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def meets(r):return r['rhat']<=1.01 and r['bulk_ess']>=400 and r['tail_ess']>=400
points=[];tasks=[]
for reg in read(out/'registry.csv'):
    tid,n=int(reg['task_index']),int(reg['n'])
    status=read(out/f'audit/task_{tid:03d}_esth_FALSE_status.csv')[0]
    task=dict(task_index=tid,n=n,status=status['status'],primary_pass=None,h_series=0,h_meeting=0,all_h_meet=None,five_h_meet=None,h_meeting_fraction=None,error_message=status['error_message'])
    if status['status']=='completed':
        d=read(out/f'audit/task_{tid:03d}_esth_FALSE_diagnostics.csv')
        h=[r for r in d if r['parameter'].startswith('h[')]
        assert len(h)==n and [int(r['parameter'][2:-1]) for r in h]==list(range(1,n+1))
        g=[]
        for r in h:
            index=int(r['parameter'][2:-1])
            q={k:float(r[k]) for k in ('mean','sd','q025','q975','rhat','bulk_ess','tail_ess','mcse_mean','mcse_sd_ratio')}
            assert all(math.isfinite(v) for v in q.values())
            g.append(dict(task_index=tid,n=n,observed_index=index,fixed_five=index in (1,n//4,n//2,3*n//4,n),primary_pass=status['primary_pass']=='TRUE',**q,meets_thresholds=meets(q)))
        passed=sum(r['meets_thresholds'] for r in g)
        task.update(primary_pass=status['primary_pass']=='TRUE',h_series=n,h_meeting=passed,all_h_meet=passed==n,five_h_meet=all(r['meets_thresholds'] for r in g if r['fixed_five']),h_meeting_fraction=passed/n)
        points.extend(g)
    tasks.append(task)
assert len(points)==29680 and len(tasks)==200
assert [r['task_index'] for r in tasks if r['status']!='completed']==[174]
assert sum(r['primary_pass'] is True for r in tasks)==22

def summarize(g):
    grouped={}
    for r in g:grouped.setdefault(r['task_index'],[]).append(r)
    fractions=[sum(r['meets_thresholds'] for r in z)/len(z) for z in grouped.values()]
    result=dict(datasets=len(grouped),h_series=len(g),series_meeting=sum(r['meets_thresholds'] for r in g),datasets_all_h_meet=sum(all(r['meets_thresholds'] for r in z) for z in grouped.values()),rhat_flagged=sum(r['rhat']>1.01 for r in g),bulk_flagged=sum(r['bulk_ess']<400 for r in g),tail_flagged=sum(r['tail_ess']<400 for r in g),dataset_fraction_median=statistics.median(fractions),dataset_fraction_min=min(fractions),dataset_fraction_max=max(fractions))
    for k in ('rhat','bulk_ess','tail_ess','mcse_sd_ratio'):
        vals=[r[k] for r in g]
        result.update({k+'_median':statistics.median(vals),k+'_min':min(vals),k+'_max':max(vals)})
    return result
summary=[];primary=[];selection=[]
for n in (40,80,160,320):
    group=[r for r in points if r['n']==n]
    for scope in ('all observed points','nested fixed five'):
        g=group if scope=='all observed points' else [r for r in group if r['fixed_five']]
        summary.append(dict(n=n,scope=scope,**summarize(g)))
        for p in (False,True):
            z=[r for r in g if r['primary_pass']==p]
            if z:primary.append(dict(n=n,scope=scope,primary_pass=p,**summarize(z)))
    ts=[r for r in tasks if r['n']==n and r['status']=='completed']
    assert not any(r['all_h_meet'] and not r['five_h_meet'] for r in ts)
    selection.append(dict(n=n,datasets=len(ts),all_observed_meet=sum(r['all_h_meet'] for r in ts),five_meet_but_other_points_flagged=sum(r['five_h_meet'] and not r['all_h_meet'] for r in ts),five_already_flagged=sum(not r['five_h_meet'] for r in ts)))
validation=read(out/'audit/validation.csv');runtime=[]
for n in (40,80,160,320):
    g=[r for r in validation if int(r['n'])==n and r['status']=='completed']
    runtime.append(dict(n=n,completed=len(g),sampling_median_seconds=statistics.median(float(r['sampling_seconds']) for r in g),postprocess_median_seconds=statistics.median(float(r['postprocess_seconds']) for r in g)))
for name,rows in [('h_by_task_point.csv',points),('by_task.csv',tasks),('h_summary.csv',summary),('h_by_primary.csv',primary),('five_vs_all.csv',selection),('runtime.csv',runtime)]:write(name,rows)
analysis=dict(completed=199,numerical_noncompletion=1,h_series=29680,primary_pass=22,summary=summary,primary=primary,selection=selection,runtime=runtime,first_draw_check_max_error=max(float(r['first_draw_check_max_error']) for r in validation if r['status']=='completed'))
analysis['independent_mean_max_error']=max(float(r['independent_mean_max_error']) for r in validation if r['status']=='completed' and r['independent_mean_max_error']!='NA')
analysis['independent_covariance_max_error']=max(float(r['independent_covariance_max_error']) for r in validation if r['status']=='completed' and r['independent_covariance_max_error']!='NA')
(out/'analysis.json').write_text(json.dumps(analysis,indent=2,allow_nan=False)+'\n')
print('PASS: 199 completed fits, 29,680 h series; nested five/all-point classification verified.')
