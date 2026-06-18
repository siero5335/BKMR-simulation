# Cross-check prompt

あなたは、BKMR 論文の最終確認を行う独立監査役として振る舞ってください。ChatGPT Pro と Claude Opus 4.8 が作成した修正文案または編集済み Word 文書を相互に比較し、どちらの出力にも残っている問題、片方だけが見落としている問題、解析結果との不整合を指摘してください。

## 比較対象

- 元文書:
  - `/Users/a_eguchi/Downloads/BKMR_sim/1_Cover_Letter.docx`
  - `/Users/a_eguchi/Downloads/BKMR_sim/2_Main_Manuscript.docx`
  - `/Users/a_eguchi/Downloads/BKMR_sim/3_Supplemental_Digital_Content.docx`
- ChatGPT Pro の修正文案または編集済みファイル
- Claude Opus 4.8 の修正文案または編集済みファイル
- 解析結果:
  - `/Users/a_eguchi/Documents/GitHub/BKMR-simulation/results/pure_bkmr_20260616_startvalues/all_task_results.csv`
  - `/Users/a_eguchi/Documents/GitHub/BKMR-simulation/results/pure_bkmr_20260616_startvalues/task_table.csv`
  - `/Users/a_eguchi/Documents/GitHub/BKMR-simulation/results/pure_bkmr_20260616_startvalues/overall_runtime.csv`
  - `/Users/a_eguchi/Documents/GitHub/BKMR-simulation/results/pure_bkmr_20260616_startvalues/parameter_diagnostics_by_task/`

## 重点監査項目

1. aBKMR、accelerated BKMR、aBKMR 由来関数、aBKMR 診断関数を使ったかのような記載が残っていないか。
2. `bkmr::SimData()`、`M = 4`、`X = 3*cos(z1) + 2*rnorm(n)` が Methods と Supplement で正しく一貫しているか。
3. `bkmr::kmbayes()` の使用、binomial/probit fit に対する明示的 `starting.values`、および「初期値指定であってモデル・事前分布・MH proposal の変更ではない」という説明が正しく入っているか。
4. Vehtari et al. (2021) に沿って、rank-normalized split R-hat と rank-normalized folded-split R-hat の最大値、bulk-ESS、tail-ESS が報告されているか。
5. R-hat < 1.01 と bulk/tail ESS >= 400 が主要基準として扱われ、R-hat < 1.05 が寛容なスクリーニングまたは感度解析として扱われているか。
6. 「431/431 fit completed, 0 errors」と「Vehtari 基準では多くが flagged」を混同していないか。
7. Results/Discussion のトーンが、収束診断の結果と整合しているか。特に「十分に収束」「安定」「問題なし」などの過剰表現が残っていないか。
8. 表、本文、補足、カバーレターで同じ数値・同じ基準が使われているか。

## 出力形式

次の形式で日本語で報告してください。

1. 重大な不整合: 採択前に必ず直すべき点。
2. 中程度の問題: 科学的トーン、Methods の精度、補足表の不足など。
3. 軽微な表現修正: 英語表現、冗長性、曖昧な言い回し。
4. ChatGPT Pro 案の良い点・弱い点。
5. Claude Opus 4.8 案の良い点・弱い点。
6. 最終的に採用すべき英文案または統合方針。

数値は必ず CSV と照合し、照合できない場合は「要確認」と明記してください。推測で数値を補わないでください。
