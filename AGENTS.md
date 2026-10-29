# AGENTS.md

Instructions for AI coding agents (Codex, Cursor, Copilot, …) are in
[CLAUDE.md](CLAUDE.md). They apply to every agent. The most important points:

- **Never upload or transmit the data** (unpublished). Never commit `data/` or rendered reports.
  The lecturer's clearance to analyse the real data with an agent covers **Claude Code** only;
  ask the lecturer before using another agent on the data.
- Use **R** / tidyverse and the helpers in `R/` (`source("R/setup.R")`).
- `A_vs_B` means `diff = log2(A) - log2(B)`. Always say which group is higher.
- Plot raw per-sample values before trusting a fold change; flag `lod_imputed` estimates.
- Cite checkable sources; keep a lab notebook in `notebook/`.
