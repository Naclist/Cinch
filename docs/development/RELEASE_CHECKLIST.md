# Cinch v0.1.0 Research Preview release checklist

Release title: **Cinch v0.1.0 — Research Preview**

| Gate | Status | Evidence |
|---|---|---|
| Naclist-owned source redistribution | PASS | `LICENSE`, `THIRD_PARTY_NOTICES.md`, `RELEASE_LICENSE_AUDIT.md` |
| Third-party source/data boundary | PASS FOR PACKAGE ARTIFACTS | no external program source vendored; SPN534 evidence excluded from sdist/wheel |
| Unified installation instructions | PASS | README explicitly checks out `integration/unified-framework` |
| Optional indexed mapper dependency | PASS | `.[mapping]` documents and installs mappy |
| Optional compiled dependency | PASS | `.[performance]` documents and installs Numba |
| Version consistency | PASS | project and CLI report `0.1.0`; frozen workflow remains `1.0.0` |
| Complete local regression | PASS | 62 tests in 3.19 s on Python 3.12; core release tests also pass on Python 3.10 |
| Tiny WGS/filter smoke | PASS | PP/PT/TP/TT 36/32/31/21; 33 final candidates |
| Conditional CLI smoke | PASS | TT differential row present; two tested hypotheses |
| sdist build and content | PASS | modules/docs/synthetic examples included; SPN534 and tests excluded |
| wheel build and content | PASS | every `cinch/*.py`, selected documentation, console entry point and license notices included |
| Twine metadata check | PASS | sdist and wheel |
| Fresh wheel installation | PASS | isolated Python 3.12 environment; mapping/performance extras; `pip check` clean |
| Fresh public GitHub clone | PASS | public HTTPS clone of `integration/unified-framework`; 62 tests, both CLI smokes, build, artifact audit, Twine check and isolated wheel install passed |
| GitHub Actions | PASS | run `36284486133` on `ae0a319`: Python 3.10, Python 3.12 and package jobs all passed |
| Production-scale acceptance | USER-DEFERRED | full SPN534, 25K-locus, and private HPC gates remain unexecuted |
| Git tag/GitHub Release | NOT CREATED | release gates pass; intentionally left for an explicit publication action |

## Scope review

[KNOWN | HIGH] Release preparation did not change a statistical definition,
state mask, permutation rule, candidate universe, threshold, mapping call, or
scientific output schema. Version reporting was separated from the frozen WGS
workflow's preserved internal version.

[KNOWN | HIGH] Passing this checklist establishes installability and Research
Preview packaging. It does not establish production readiness or biological
causality.
