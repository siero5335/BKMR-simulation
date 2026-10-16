"""Isolated rerun entry point. Python 3.9+; standard library only."""
import argparse, concurrent.futures, csv, json, os, pathlib, shutil, subprocess, sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def run_r(work, script, *args, log=None):
    rscript=shutil.which(os.environ.get('RSCRIPT','Rscript'))
    if not rscript:raise RuntimeError('Rscript not found. Set RSCRIPT to its executable path.')
    env=os.environ.copy()
    # Preserve the historical threading environment; record it instead of changing it.
    if log:
        with log.open('w') as f:
            subprocess.run([rscript,str(script),*map(str,args)],cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    else:subprocess.run([rscript,str(script),*map(str,args)],cwd=work,env=env,check=True)
def prepare(work):
    if work.exists():raise RuntimeError('prepare requires a new directory; existing files will not be replaced.')
    work.mkdir(parents=True)
    shutil.copytree(ROOT/'code/runtime',work,dirs_exist_ok=True)
    shutil.copytree(ROOT/'reference',work/'reference')
    shutil.copytree(ROOT/'environment',work/'environment')
    (work/'code').mkdir(exist_ok=True)
    for n in ('validate_inputs.R','prepare_initialization.R','check_environment.R'):
        shutil.copy2(ROOT/'code'/n,work/'code'/n)
    res=work/'Main/Code_Data/results';res.mkdir(parents=True)
    for n in ('task_table.csv','all_task_results.csv'):shutil.copy2(ROOT/'coauthor_additions'/n,res/n)
    d=res/'parameter_diagnostics_by_task';d.mkdir()
    groups={}
    for r in read(ROOT/'reference/original_parameter_diagnostics.csv'):
        tid=int(r.pop('task_index'));groups.setdefault(tid,[]).append(r)
    for tid,rows in groups.items():
        with (d/f'task_{tid:03d}__reference.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for stem in ('revision_20260911','h_false_allpoints_20260912','initialization_extension_20260914'):
        for sub in ('audit','logs','checkpoints','figures','failure_case','chains','inputs'):
            (work/stem/sub).mkdir(parents=True,exist_ok=True)
    for n in ('registry.csv','cases.csv','plan.json'):
        shutil.copy2(ROOT/'initialization_extension'/n,work/'initialization_extension_20260914'/n)
    # Figure regeneration can run without MCMC using copies of published CSVs.
    for n in ('fit_summaries.csv','parameter_diagnostics.csv','beta_by_chain.csv','chain_status.csv'):
        shutil.copy2(ROOT/'initialization_extension'/n,work/'initialization_extension_20260914'/n)
    rows=[]
    for t in read(res/'task_table.csv'):
        tid=int(t['task_index'])
        if 51<=tid<=250:
            rows.append(dict(task_index=tid,n=t['n'],reused='False',source='h_false_allpoints_20260912',
                runner='h_false_allpoints_20260912/code/run_task_v2.R',checkpoint=f'h_false_allpoints_20260912/checkpoints/task_{tid:03d}_esth_FALSE.rds'))
    with (work/'h_false_allpoints_20260912/registry.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (work/'WORKSPACE.json').write_text(json.dumps({'prepared_from':str(ROOT),'kind':'BKMR Rev1 isolated rerun','published_csv_copies':'Initialization summary CSVs initially copied for figure-only regeneration; init-summary replaces them only after validating all 144 chain files.'},indent=2)+'\n')
    print('Prepared',work)
def batch(work,stage,ids,workers):
    script={'beta':'revision_20260911/code/02_reproduce_task.R','h':'h_false_allpoints_20260912/code/run_task_v2.R','init':'initialization_extension_20260914/code/run_chain.R'}[stage]
    logdir=work/'run_logs';logdir.mkdir(exist_ok=True)
    # Stop new dispatch on failure; already running jobs finish.
    todo=iter(ids);failures=[]
    def one(i):run_r(work,script,i,log=logdir/f'{stage}_{i}.log');return i
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        active={}
        for _ in range(min(workers,len(ids))):
            i=next(todo);active[pool.submit(one,i)]=i
        while active:
            done,_=concurrent.futures.wait(active,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in done:
                i=active.pop(f)
                try:f.result();print('Completed',stage,i,flush=True)
                except Exception as exc:failures.append((i,str(exc)));print('FAILED',stage,i,flush=True)
            if not failures:
                while len(active)<workers:
                    i=next(todo,None)
                    if i is None:break
                    active[pool.submit(one,i)]=i
    if failures:raise RuntimeError(f'{failures}; inspect {logdir}')
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['prepare','check','audit','beta','beta-summary','h','h-summary','init-inputs','init','init-summary','comparisons','figures','failure-trace','baseline-summary','initial-values'])
    p.add_argument('--workdir',required=True,type=pathlib.Path)
    select=p.add_mutually_exclusive_group();select.add_argument('--task',type=int);select.add_argument('--job');select.add_argument('--all',action='store_true');select.add_argument('--init-cases',action='store_true')
    p.add_argument('--workers',type=int,default=1)
    a=p.parse_args();work=a.workdir.expanduser().resolve()
    if a.workers<1:p.error('--workers must be positive')
    if a.stage=='prepare':prepare(work);return
    if not (work/'WORKSPACE.json').is_file():p.error('Run prepare first with a new directory.')
    if a.stage=='check':
        run_r(work,'code/check_environment.R','environment/package_versions.csv');run_r(work,'code/validate_inputs.R');return
    if a.stage=='init-inputs':
        run_r(work,'code/prepare_initialization.R',*([a.task] if a.task is not None else []));return
    if a.stage in ('beta','h','init'):
        if a.stage=='init':
            jobs=[r['job_id'] for r in read(work/'initialization_extension_20260914/registry.csv')]
            if a.all:ids=jobs
            elif a.job in jobs:ids=[a.job]
            else:p.error('init requires --all or a valid --job from registry.csv')
        else:
            if a.all:ids=list(range(51,251))
            elif a.init_cases:ids=[int(r['task_index']) for r in read(work/'initialization_extension_20260914/cases.csv')]
            elif a.task is not None and 51<=a.task<=250:ids=[a.task]
            else:p.error('beta/h requires --task 51..250, --init-cases, or --all')
            if a.stage=='h':
                for i in ids:
                    if not (work/f'revision_20260911/checkpoints/task_{i:03d}.rds').exists():p.error(f'Run beta for task {i} first.')
        batch(work,a.stage,ids,a.workers);return
    scripts={
        'audit':['revision_20260911/code/01_audit_existing.R'],
        'beta-summary':['revision_20260911/code/03_summarize_beta.R'],
        'h-summary':['h_false_allpoints_20260912/code/validate_full.R','h_false_allpoints_20260912/code/summarize.py'],
        'init-inputs':['code/prepare_initialization.R'],
        'init-summary':['initialization_extension_20260914/code/summarize.R','initialization_extension_20260914/code/summarize_comparisons.py'],
        'comparisons':['initialization_extension_20260914/code/summarize_comparisons.py'],
        'figures':['initialization_extension_20260914/code/make_figures.R'],
        'failure-trace':['revision_20260911/code/04_trace_failure.R'],
        'baseline-summary':['initialization_extension_20260914/code/baseline_beta.R'],
        'initial-values':['initialization_extension_20260914/code/export_initial_values.R']}
    for script in scripts[a.stage]:
        if script.endswith('.py'):subprocess.run([sys.executable,script],cwd=work,check=True)
        else:run_r(work,script)
if __name__=='__main__':main()
