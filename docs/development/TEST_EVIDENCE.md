# Test evidence

| Date | Source/commit | Environment | Command | Result | Scope |
|---|---|---|---|---|---|
| 2026-09-26 | CINCH `98c41ef` | Python 3.12.14 local venv | `pytest -q` | 12 passed in 5.88 s | core semantics and frozen release assets |
| 2026-09-26 | CINCH-dev2 `ec1eaa5` | Python 3.12.14 dev2 venv | `pytest -q tests` | 7 passed | weighted MI/EpiDis/BH and frozen hashes |
| 2026-09-26 | CINCH `98c41ef` | system Python 3.9.6 | editable install | failed before tests | environment below declared `>=3.10`; not a code failure |

## Evidence limitations

[KNOWN | HIGH] Existing tests do not validate legacy uberBlast installation, true mapping parallelism, resumability, indels/truncations/competing loci, unified numerical equivalence, 25K-locus memory, or end-to-end integration.
