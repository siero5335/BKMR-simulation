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

methods=('To examine the practical implications of initialization and longer computation, we selected three completed est.h = FALSE datasets at each of n = 40, 80, 160, and 320, with variable selection disabled. '
 'Selection used baseline diagnostics and was fixed before the extension runs: the previously used first-task anchor, the lowest-index primary-pass case, and the remaining primary-non-pass case with the largest beta R-hat. '
 'Because there was no primary-pass case at n = 40, we used the middle-ranked beta R-hat case among the remaining datasets instead. These 12 cases were illustrative rather than a representative simulation sample. '
 'For each dataset, we compared the original four starting states with two newly generated sets of four starts from the same initialization distributions. The four sampler seeds were held fixed across initialization sets; data, priors, proposals, and the FALSE sampling algorithm were unchanged. '
 'Each of the resulting 144 chains was run toward 100,000 iterations. We evaluated prefixes ending at 10,000, 20,000, 40,000, and 100,000 iterations, primarily excluding the first 5,000 iterations. '
 'A sensitivity analysis excluded the first half of each prefix. The original initialization runs exactly reproduced iterations 5001–10000 of all six parameter series from the previous original-condition fits in all 48 chains. '
 'We evaluated beta means, 95% equal-tailed intervals, errors relative to the generating value 0.5, interval inclusion of zero and the true value, and MCSEs, alongside the unchanged six-parameter R-hat/bulk-ESS/tail-ESS criterion. '
 'All six-parameter trajectories were retained. Numerical failures were recorded separately; a prefix was evaluated only if all four chains had completed it, excluding any partially executed failing iteration. Later failures remained documented when shorter prefixes were evaluable.')
results=(f'Of the 144 chains, {144-len(errors)} reached 100,000 iterations and {len(errors)} returned sampler errors. '
 f'At 10,000 iterations, {pass10}/{len(ok10)} evaluable dataset–initialization configurations met the six-parameter primary criterion; at 100,000 iterations, {pass100}/{len(full_ok)} did so. '
 f'Among the {len(pair)} configurations evaluable at both endpoints, {trans["nonpass_to_pass"]} changed from non-pass to pass, {trans["pass_to_nonpass"]} from pass to non-pass, {trans["both_pass"]} passed at both, and {trans["both_nonpass"]} passed at neither. '
 f'Among the {len(common_ids)} datasets with all three initialization sets evaluable at both endpoints, the median range of beta posterior means across initialization sets was {fmt(median(before))} at 10,000 iterations and {fmt(median(after))} at 100,000 iterations. '
 f'This range decreased in {improved} datasets and increased in {increased}. '
 f'At 100,000 iterations, the fixed-5,000 and half-chain warmup summaries differed in primary status for {notes["warmup_fixed_only"]+notes["warmup_half_only"]} configurations ({notes["warmup_fixed_only"]} passed only with fixed warmup; {notes["warmup_half_only"]} only with half-chain warmup).')
limits=('These comparisons assess initialization sensitivity and the practical effect of staged computation in selected cases. Repeated initialization sets use the same datasets and paired sampler seeds; they are not independent simulation replicates or independent sets of chains to be pooled for ESS estimation. '
 'The prefixes and warmup sensitivity summaries also reuse the same trajectories. We therefore report descriptive changes, not population convergence rates, causal effects of diagnostic failure, or a new coverage study. '
 'Intervals and tail probabilities from diagnostic-non-pass fits are exploratory summaries of computed draws. Small MCSEs do not establish adequate posterior exploration. '
 'Increased stability across starts does not by itself establish accuracy relative to the generating value. Persistent difficulties at 100,000 iterations document an unresolved limitation under the tested model, sampler, and computational budget; weak identification is a possible explanation alongside algorithmic and numerical limitations, not a conclusion established by non-pass alone. '
 'No new h, exposure-response curve, or mixture-effect contrast assessment was performed in this extension; the observed-point h results remain a separate supplementary analysis.')
table=['### Diagnostic outcomes by iteration budget','',
 '| n | Iterations per chain | Evaluable configurations / planned | Six-parameter primary pass / evaluable | Beta-only diagnostic pass / evaluable | Median beta MCSE/SD |',
 '|---|---:|---:|---:|---:|---:|']
for r in summary:table.append(f"| {r['n']} | {r['budget']:,} | {r['evaluable']}/{r['configurations']} | {r['primary_pass']}/{r['evaluable']} | {r['beta_pass']}/{r['evaluable']} | {fmt(100*r['beta_mcse_sd_median'],2)}% |")
table+=['','Each n includes three selected datasets with three initialization sets each. These nine configurations are dependent, descriptive comparisons. Warmup is fixed at 5,000.','',
 '### Change in initialization sensitivity within datasets','',
 '| n | Task | Selection role | All 3 sets evaluable at 100k | Beta mean range across starts: 10k | Beta mean range across starts: 100k | Primary-pass sets at 100k / evaluable |',
 '|---|---:|---|---|---:|---:|---:|']
for c in cases:
    tid=int(c['task_index']);x,y=ilook[tid,10000],ilook[tid,100000]
    table.append(f"| {c['n']} | {tid} | {c['selection_role']} | {y['all_three_evaluable']} | {fmt(x['beta_mean_range'])} | {fmt(y['beta_mean_range'])} | {y['primary_pass_sets']}/{y['evaluable_sets']} |")
table+=['','Ranges are not reported when any of the three initialization sets is unavailable. Figures show beta means and 95% intervals, maximum six-parameter R-hat, and minimum bulk/tail ESS; the generating beta is 0.5. Full configuration- and chain-specific summaries are supplied as CSV files.']
supp='\n\n'.join([methods,results,limits,'\n'.join(table)])+'\n'
(out/'SUPPLEMENT_ADDITION.md').write_text('# Initialization and chain-length sensitivity\n\n'+supp)
ja=['# 共著者提案を反映した追加実験：初期値と反復数','',
 '初期値を変えたときに推論結果が変わるか、反復数を増やすと改善するかを実際に調べました。n=40/80/160/320から各3件、計12件を選び、各データセットで初期値を3組×4鎖、最大10万反復まで計算しました。元のFALSEサンプラー、データ、事前分布、proposalは保持しています。','',
 f'144本の鎖のうち{144-len(errors)}本が10万反復に到達し、{len(errors)}本でsampler errorが記録されました。元の初期値の48本は、最初の1万反復の事後部分が従来の元6パラメータdrawsと完全一致しました。','',
 '## 反復数を増やすとどうなったか','',
 f'1万反復では評価可能な{len(ok10)}条件中{pass10}条件、10万反復では{len(full_ok)}条件中{pass100}条件が元6パラメータのprimary基準を満たしました。分母が異なる場合は数値未完了によるものです。','',
 f'両時点を比較できた{len(pair)}条件では、不通過→通過が{trans["nonpass_to_pass"]}、通過→不通過が{trans["pass_to_nonpass"]}、両方通過が{trans["both_pass"]}、両方不通過が{trans["both_nonpass"]}でした。','',
 '| n | 反復数/鎖 | 評価可能/予定条件 | primary通過/評価可能 | β MCSE/SD中央値 |',
 '|---|---:|---:|---:|---:|']
for r in summary:ja.append(f"| {r['n']} | {r['budget']:,} | {r['evaluable']}/9 | {r['primary_pass']}/{r['evaluable']} | {fmt(100*r['beta_mcse_sd_median'],2)}% |")
ja+=['','各nの9条件は3データセット×初期値3組です。独立した9データセットではなく、通過率をBKMR一般の成功率と解釈できません。primary比較のwarmupは5,000に固定しました。','',
 '## 初期値で推論結果はどれだけ変わったか','',
 f'両時点ですべての初期値組を比較できた{len(common_ids)}データセットでは、初期値組間のβ推定平均の最大値−最小値の中央値は、1万反復で{fmt(median(before))}、10万反復で{fmt(median(after))}でした。範囲が小さくなったのは{improved}件、大きくなったのは{increased}件です。',
 'βの推定平均・95%区間と真値0.5を各ケースで比較した図を添付しています。真値との差と、初期値間の再現性は区別して確認します。基準不通過の区間・確率は、得られたdrawsの探索的要約です。','',
 '## warmupの取り方も確認','',
 f'10万反復時に、最初の5,000反復を除く集計と最初の50,000反復を除く集計でprimaryが異なったのは{notes["warmup_fixed_only"]+notes["warmup_half_only"]}条件でした。固定warmupだけが通過したのは{notes["warmup_fixed_only"]}、後半だけの集計が通過したのは{notes["warmup_half_only"]}です。都合のよい集計に置き換えず、両方を示します。','',
 '## 共著者の提案をどう反映するか','',
 '初期値を変えた推論結果の比較と、反復数を段階的に増やした実験をSupplementに追加し、AE-1とComment 6の回答に反映しました。既存のβのbias・coverage解析は、真値に対する推定性能を直接調べる中心的な解析として残しています。',
 '実務上は、初期値間の違いと反復数延長後の診断・MCSE・推定値の変化を確認し、改善があれば延長する、改善が乏しければ計算方法や識別の弱さを調査する、という判断手順にします。10万反復という値は今回の検討範囲であり、普遍的な必要反復数ではありません。',
 '長く回しても問題が残る場合は、検討した条件・計算予算では十分な診断上の根拠が得られなかったと報告します。データの情報量不足・識別の弱さは原因候補として考察し、数値実装やサンプラーの影響と区別します。','',
 '## 選定と保存','',
 '既存の先行ケース、primary通過例、βのR-hatが特に高い不通過例を、今回の延長結果を見る前に規則で選びました。n=40にはprimary通過例がないため、中央順位のβ R-hat例で代用しています。選定された12件の例示であり、標本サイズ間や全データセットの比較を代表するものではありません。',
 '初期値3組は元と同じ生成範囲を使い、各鎖のsampler seedを組間でそろえています。組間・反復数間の結果は依存するため、独立した反復として扱っていません。今回の全144本の元6パラメータのraw drawsはローカルに保存しています。hは今回再計算していません。']
(out/'RESULTS_JA.md').write_text('\n'.join(ja)+'\n')

base_dir=Path('deliverables/BKMR_manuscript_FALSE_h_allpoints_20260913')
letter=(base_dir/'02_Response_letter.md').read_text()
addition=('We additionally examined sensitivity to initialization and chain length in 12 selected datasets, using three sets of four initial states and prefixes of 10,000, 20,000, 40,000, and 100,000 iterations under the unchanged FALSE model and sampler. '
 +results+' These illustrative comparisons supplement the truth-based beta analysis above and are detailed in [Supplemental section: initialization and chain-length sensitivity].\n\n')
marker='### AE-2. Clarify that h was neither sampled nor monitored'
assert letter.count(marker)==1;letter=letter.replace(marker,addition+marker)
marker2='The template is a reporting aid, not a new convergence criterion, and no acceptance-rate threshold is proposed.'
assert letter.count(marker2)==1
letter=letter.replace(marker2,
 'We also evaluated staged iteration extension and initialization sensitivity empirically in the selected cases, rather than recommending additional iterations without examining their practical effect. '
 +f'Among {len(pair)} configurations evaluable at both 10,000 and 100,000 iterations, {trans["nonpass_to_pass"]} changed from primary non-pass to pass and {trans["pass_to_nonpass"]} changed in the reverse direction. '
 'The workflow now explicitly recommends comparing inferential summaries across dispersed starts, assessing whether added iterations improve diagnostics and MCSE, and reporting unresolved problems when the tested budget is exhausted. '
 'Weak identification is considered a possible contributor to persistent difficulties, alongside computational and numerical limitations; non-pass after a long run is not used alone to diagnose insufficient information in the data.\n\n'+marker2)
ms=(base_dir/'03_Manuscript_and_Supplement.md').read_text()
short=('## M13. Initialization and iteration sensitivity — Methods\n\n'
 'In 12 cases selected using baseline diagnostics, we compared three sets of four initial states under the unchanged FALSE model and sampler and evaluated 10,000-, 20,000-, 40,000-, and 100,000-iteration prefixes. We assessed beta inferential summaries and the original six-parameter diagnostics, with fixed-warmup and half-chain sensitivity summaries (Supplementary Information).\n\n'
 '## M14. Initialization and iteration sensitivity — Results\n\n'
 +f'Among {len(pair)} dataset–initialization configurations evaluable at both 10,000 and 100,000 iterations, {trans["nonpass_to_pass"]} changed from primary non-pass to pass and {trans["pass_to_nonpass"]} changed from pass to non-pass. '
 +f'The median range of beta posterior means across initializations changed from {fmt(median(before))} to {fmt(median(after))} among {len(common_ids)} datasets with all three initialization sets evaluable at both endpoints. '
 'Numerical non-completion and persistent diagnostic limitations were retained explicitly (Supplementary Information).\n\n')
assert ms.count('## S1.')==1;ms=ms.replace('## S1.',short+'## S1.',1)
workflow=('Analysts should compare chains initialized at dispersed starting points and assess both diagnostics and the inferential summaries used for interpretation. '
 'If mixing is satisfactory but ESS or MCSE indicates insufficient precision, increase iterations and reassess. If chain disagreement persists, extend in stages while checking whether diagnostics, MCSE, and estimates improve; lack of improvement warrants investigation of parameterization, sampling efficiency, numerical implementation, and possible weak identification. '
 'The selected-case extension to 100,000 iterations illustrates this evaluation and does not establish a universal iteration requirement. '
 'Changes to priors alter the posterior target and require scientific justification and separate sensitivity analysis. Retain all chain outcomes and numerical failures. '
 'When adequate diagnostic evidence cannot be obtained within the tested budget, report the unresolved inferential limitation; the result alone does not identify insufficient data information as its cause.')
ms,n=re.subn(r'Analysts should examine multiple chains initialized at dispersed points.*?When adequate diagnostic evidence cannot be obtained, report the unresolved limitation\.',lambda _:workflow,ms,count=1,flags=re.S);assert n==1
ms+='\n## S8. Initialization and chain-length sensitivity\n\n'+supp
(out/'drafts/02_Response_letter_filled.md').write_text(letter)
(out/'drafts/03_Manuscript_and_Supplement_filled.md').write_text(ms)
subprocess.run(['/usr/local/bin/Rscript',str(out/'code/make_figures.R')],check=True)

stage=Path('deliverables/BKMR_manuscript_initialization_extension_20260914')
if not stage.exists():shutil.copytree(base_dir,stage)
shutil.copy2(out/'drafts/02_Response_letter_filled.md',stage/'02_Response_letter.md')
shutil.copy2(out/'drafts/03_Manuscript_and_Supplement_filled.md',stage/'03_Manuscript_and_Supplement.md')
extra_dir=stage/'additional_initialization_extension';extra_dir.mkdir(exist_ok=True)
for name in ('RESULTS_JA.md','SUPPLEMENT_ADDITION.md','PLAN_JA.md','PROVENANCE_JA.md','cases.csv','initial_values.csv','budget_summary.csv','initialization_sensitivity.csv','fit_summaries.csv','parameter_diagnostics.csv','beta_by_chain.csv','chain_status.csv','baseline_beta_by_chain.csv','baseline_beta_by_task.csv','plan.json'):
    shutil.copy2(out/name,extra_dir/name)
for name in ('paired_extension.csv','warmup_sensitivity.csv'):
    if (out/name).exists():shutil.copy2(out/name,extra_dir/name)
shutil.copytree(out/'figures',extra_dir/'figures',dirs_exist_ok=True)
shutil.copytree(out/'code',extra_dir/'code',dirs_exist_ok=True)
(stage/'00_README_JA.md').write_text('''# 共著者提案を反映した原稿作成用の最新版

最初に `additional_initialization_extension/RESULTS_JA.md` をご覧ください。12データセットで初期値3組×4鎖、最大10万反復を検討した追加結果です。AE-1、Comment 6、本文追記案、Supplement S8を更新しました。既存のbetaの推定性能、FALSEの全観測点h、TRUEとの比較、数値エラー解析は保持しています。

新しい実行・集計コードと、全条件の診断CSV・鎖別beta要約・図を収録しています。コードは実行時の作業フォルダ構成を参照するため、元環境の外でそのまま動作する持ち運び用一式ではありません。選定例の結果を全体の成功率や情報量不足の証明として解釈しないでください。hは今回再計算していません。

原稿用ZIPには今回の大容量のraw drawsを含めていません。元6パラメータの全144本の軌跡は作業フォルダ `initialization_extension_20260914/chains/` に保存しています。別の検証用ZIPには、全12件の入力と、事前に選んだtask 201の元の初期値組・alternative_1の計8本のraw trajectoriesを収録しています。

DOCXへの直接挿入、Paperpile、表番号・ページ行番号の確定は未実施です。レターの完了形は原稿へ実際に挿入した後に使用してください。
''')

verify=Path('deliverables/BKMR_initialization_extension_verification_20260914')
verify.mkdir(exist_ok=True)
for directory in ('inputs','code'):
    shutil.copytree(out/directory,verify/directory,dirs_exist_ok=True)
(verify/'raw_chains').mkdir(exist_ok=True)
for set_name in ('original','alternative_1'):
    for chain in range(1,5):
        name=f'task_201_{set_name}_chain_{chain}.rds';shutil.copy2(out/'chains'/name,verify/'raw_chains'/name)
for name in ('PROVENANCE_JA.md','cases.csv','initial_values.csv','registry.csv','plan.json','fit_summaries.csv','parameter_diagnostics.csv','chain_status.csv'):
    shutil.copy2(out/name,verify/name)
shutil.copytree(out/'audit',verify/'audit',dirs_exist_ok=True)
(verify/'README_JA.md').write_text('''# 初期値・反復数実験の検証用追補

全12件の入力・初期値・seed・元の6パラメータの参照drawsと、実際の実行・集計コード、全条件の診断CSVを収録しました。raw_chainsにはtask 201（n=40）のoriginalとalternative_1、各4鎖、計8本のraw trajectoriesを収録しています。初期反復・warmupも含むため、診断時には該当する反復だけを選んでください。

RDSのdrawsは行が反復1からcompleted_iterations、列がbeta[1], lambda, r[z1]〜r[z4]です。短い軌跡はsampler_errorの可能性があるためstatusとerror_iterationを確認してください。10k/20k/40k/100k時点の主比較は行5001:budget、感度分析はfloor(budget/2)+1:budgetを使用します。各初期値組内の4鎖で診断し、異なる初期値組を独立な8鎖としてプールしないでください。

独立したraw drawsからの再集計は、展開したフォルダ内で `Rscript code/recompute_raw_diagnostics.R . NEW_AUDIT_DIRECTORY` を実行します。posterior 1.7.0を使い、今回収録した2初期値組について、評価可能な全反復数・warmupの9指標をCSVと照合します。サンプラーは再実行しません。\n\nこの部分集合のraw drawsで全144本を独立検証したとは言えません。全raw drawsはローカルに保持しています。入力・コードには当時の作業フォルダとR/package fingerprintへの参照があるため、別環境でのそのままの実行は保証しません。原稿用の成果物は別のZIPです。
''')
with tempfile.TemporaryDirectory(prefix='bkmr_extension_audit_') as temporary:
    audit_output=Path(temporary)/'new_audit'
    subprocess.run(['/usr/local/bin/Rscript',str(out/'code/recompute_raw_diagnostics.R'),str(verify),str(audit_output)],check=True)
    shutil.copytree(audit_output,verify/'audit/raw_recomputation',dirs_exist_ok=True)
    shutil.copytree(audit_output,out/'audit/raw_recomputation_task201',dirs_exist_ok=True)
def package(folder):
    archive=folder.with_suffix('.zip');temp=archive.with_suffix('.zip.tmp')
    files=sorted(p for p in folder.rglob('*') if p.is_file())
    with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in files:z.write(p,arcname=(Path(folder.name)/p.relative_to(folder)).as_posix())
    with zipfile.ZipFile(temp) as z:
        assert z.testzip() is None
        for p in files:assert z.read((Path(folder.name)/p.relative_to(folder)).as_posix())==p.read_bytes()
    temp.replace(archive);sha=hashlib.sha256(archive.read_bytes()).hexdigest()
    archive.with_suffix('.zip.sha256').write_text(sha+'  '+archive.name+'\n')
    return dict(archive=str(archive.resolve()),bytes=archive.stat().st_size,files=len(files),sha256=sha)
delivery=dict(manuscript=package(stage),verification=package(verify))
(out/'delivery.json').write_text(json.dumps(delivery,indent=2)+'\n')
state=json.loads((out/'batch_status.json').read_text());elapsed=time.time()-state['started_unix']
(out/'COMPLETION.md').write_text(f'# 完了\n\n144 chain jobs; elapsed through packaging {elapsed/3600:.2f} hours. Completed chains: {144-len(errors)}; sampler errors: {len(errors)}. All 48 original 10k prefixes verified.\n')
print(json.dumps(delivery,indent=2))
