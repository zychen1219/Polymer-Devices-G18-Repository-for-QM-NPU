# Freeze record

Written when the pipeline was frozen, before submission.

| Item | Value |
|---|---|
| Pipeline commit (contains `src/`, `tests/`, `tools/`, `data/`) | `2008bef947f8bbf792a317c3c180ab9d2e95b80f` |
| Repository HEAD when this record was written | `2daec63` (the pin commit, which follows the pipeline commit) |
| Pinned in | `results/results.json` (`pipeline_commit`), `blind/pipeline_commit.txt` |
| Freeze date | 2026-11-17 |

## How the two hashes relate

`reproduce.py` and `run_blind.py` read `git rev-parse HEAD` when they run and write that
value into their outputs. A run made from `2008bef` therefore records `2008bef`. Writing
that value into the repository then requires a commit, which moves HEAD forward, so the
recorded hash is always an ancestor of HEAD rather than HEAD itself. Re-running the
pipeline at a later commit does not change any number or figure; it only replaces the
recorded hash with the newer one. The value above is the commit that produced the numbers
in `results/results.json`. Later documentation commits move HEAD further forward; the pipeline
hash above stays valid because the code it names does not change.

## Commands that produced the frozen outputs

```bash
python reproduce.py
python run_blind.py --data data/raw/blind_G18 --out blind
python -m unittest discover -s tests
```

Requirements: Python 3.12 with numpy and matplotlib (`requirements.txt`). The vendored
wheel directory and the matplotlib font cache are local conveniences and are excluded by
`.gitignore`; a clean clone installs nothing but numpy and matplotlib.
