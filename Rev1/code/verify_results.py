"""Verify distributed CSV summaries without MCMC. Python stdlib + base R."""
from pathlib import Path
import csv, json, math, subprocess, sys, tempfile, shutil, os
ROOT=Path(__file__).resolve().parents[1]
def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def equal(a,b):
    if a==b:return True
    try:
        x,y=float(a),float(b)
        return (math.isnan(x) and math.isnan(y)) or math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-12)
    except (ValueError,TypeError):return False
def compare(a,b):
    aa,bb=rows(a),rows(b)
    assert len(aa)==len(bb),(a,b,'row count')
    for i,(x,y) in enumerate(zip(aa,bb)):
        assert x.keys()==y.keys(),(a,b,'columns')
        assert all(equal(x[k],y[k]) for k in x),(a,b,i)
R=shutil.which(os.environ.get('RSCRIPT','Rscript'))
if not R:raise SystemExit('Base R is required for the independent h reaggregation.')
with tempfile.TemporaryDirectory(prefix='bkmr-rev1-verify-') as temp:
    dest=Path(temp)
    subprocess.run([sys.executable,str(ROOT/'coauthor_additions/reproduce_summaries.py'),'--output-dir',str(dest/'coauthor')],check=True,stdout=subprocess.PIPE,text=True)
    assert json.loads((dest/'coauthor/numeric_checks.json').read_text())==json.loads((ROOT/'coauthor_additions/numeric_checks.json').read_text())
    for n in ('reference_comparisons.csv','chain_ranges.csv'):compare(dest/'coauthor'/n,ROOT/'coauthor_additions'/n)
    # Run the extracted historical comparisons against copied summary inputs.
    init=dest/'initialization_extension_20260914';init.mkdir()
    for n in ('fit_summaries.csv','chain_status.csv','cases.csv'):shutil.copy2(ROOT/'initialization_extension'/n,init/n)
    subprocess.run([sys.executable,str(ROOT/'code/runtime/initialization_extension_20260914/code/summarize_comparisons.py')],cwd=dest,check=True)
    for n in ('initialization_sensitivity.csv','budget_summary.csv','paired_extension.csv','warmup_sensitivity.csv'):compare(init/n,ROOT/'initialization_extension'/n)
    subprocess.run([R,str(ROOT/'code/reaggregate_h.R'),str(dest/'h')],check=True)
    # Independently rederive six-parameter flags from 430 task-level diagnostic sets.
    grouped={}
    for r in rows(ROOT/'reference/original_parameter_diagnostics.csv'):
        grouped.setdefault(int(r['task_index']),[]).append(r)
    original={int(r['task_index']):r for r in rows(ROOT/'coauthor_additions/all_task_results.csv')}
    assert len(original)==431 and len(grouped)==430
    passes=0
    for tid,g in grouped.items():
        assert len(g)==6
        ok=all(math.isfinite(float(r[k])) for r in g for k in ('rhat','bulk_ess','tail_ess')) and all(float(r['rhat'])<=1.01 and float(r['bulk_ess'])>=400 and float(r['tail_ess'])>=400 for r in g)
        assert ok==(original[tid]['converged_hard'].lower()=='true')
        passes+=ok
    assert passes==30
print('PASS: coauthor summaries, four initialization comparison tables, all 29,680 h point records, and original 430 primary classifications. No MCMC performed; distributed results unchanged.')
