# Document revision prompts

These prompts are intended for independent review with ChatGPT Pro and Claude
Opus 4.8.

Input Word files:

- `/Users/a_eguchi/Downloads/BKMR_sim/1_Cover_Letter.docx`
- `/Users/a_eguchi/Downloads/BKMR_sim/2_Main_Manuscript.docx`
- `/Users/a_eguchi/Downloads/BKMR_sim/3_Supplemental_Digital_Content.docx`

Pure BKMR result folder:

- `/Users/a_eguchi/Documents/GitHub/BKMR-simulation/results/pure_bkmr_20260616_startvalues/`

Use:

1. Run `01_prompt_for_chatgpt_pro.md` in ChatGPT Pro.
2. Run `02_prompt_for_claude_opus_4_8.md` in Claude Opus 4.8.
3. After both return revised text or edited documents, use
   `03_cross_check_prompt.md` to compare the two outputs.

The prompts intentionally require the models to distinguish fit completion from
MCMC convergence diagnostics and to avoid unsupported claims of convergence.
