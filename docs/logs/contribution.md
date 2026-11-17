# Contribution statement

Shares were agreed before submission and match the authorship of the commits in this repository.
Each member owns a set of files, reviewed the code or text in that set, and can explain it.

| Member | Student ID | Share | Files owned (commit author) | Work done |
|---|---|---|---|---|
| Chen Zhengyi | 2023304324 | 34% | `src/common.py`, `src/l1.py`, `src/l3.py`, `reproduce.py`, report sections 2 and 4 | Fixed the L1 model definition and the working-point convention; built the L3 reconstruction, the sensitivity sweep and the noise floor; integrated the pipeline so that one command regenerates every figure and `results.json`; wrote the L1 and L3 text |
| Jiang Tianyu | 2023304369 | 30% | `src/l2.py`, `run_blind.py`, report section 3 | Built the cone/gap model pair and the nearest-point fit; set the separate normalisation and the 15 per cent selection rule; ran and documented the 200-draw resampling; wrote the L2 text and ran the blind pipeline |
| Qiao Qinpeng | 2023304299 | 19% | `tests/`, `data/raw/`, `figures/` | Assembled the data tree and the manifest checks; wrote the assertion suite that caught the sign and zero-mean errors; produced and checked the seven figures; re-ran every stage after each fix |
| Wang Ziyue | 2023304298 | 17% | `docs/logs/`, `docs/references.bib`, report section 5 and the reference list | Kept the AI usage log and the error ledger; verified the citations and DOIs; checked the word count, the figure order and the cross-references; drafted and revised the Transfer section |

The percentages reflect both the amount of work and its weight in the mark scheme: the two
computational sections carry 35 of the 100 group marks and the written argument carries 10, so a
larger share of the model-building work is a larger share of the submission.

Working practice: the group met at the start and end of each stage. Models and acceptance criteria
were set before any code was written; the assistant drafts described in `ai-usage-log.md` were
accepted only after the gates in `tests/` passed. Disagreements about the sign convention and about
the choice of eta were settled by re-running the pipeline, not by argument, and both are recorded
in `error-ledger.md`.
