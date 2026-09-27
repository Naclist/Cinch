# Cinch v0.1.0 — Research Preview

## Release status

Cinch v0.1.0 is the first unified Research Preview. It is intended for
reproducible research evaluation and method development. It is not a
production-ready clinical, diagnostic, epidemiological, or regulatory tool.

The immutable release candidate is prepared on
`integration/unified-framework`. A Git tag and GitHub Release must be created
only after the final public-clone, wheel, smoke-test, and GitHub Actions checks
pass.

## Validation classes

### Validated functionality

- The frozen WGS and filter implementations remain available without changing
  their scientific definitions.
- Weighted binary MI, EpiDis, and BH match the preserved Cinch-dev2 kernels on
  controlled equivalence tests.
- PROFILE_V2 conversion, blockwise PP/PT/TP/TT association, sparse distance
  handling, resumability, staged PP SHC/BH/ARACNE filtering, and report
  generation pass the repository regression suite.
- The indexed mapper has controlled tests for cache identity, parallel genome
  execution, allele calls, and resume behavior. It is not claimed to reproduce
  historical uberBlast allele nomenclature.

### Controlled-method validation

- Population-conditioned and two-background differential association are
  implemented for PP, PT, TP, and multiallelic TT.
- Cases A–J cover population-only association, habitat-specific coexistence and
  exclusion, a null contrast, perfect habitat-lineage confounding,
  multiallelic states, PT/TP conditioning, missing/type ambiguity, unbalanced
  support, and preserved frozen Diff-GWES behavior.
- These controlled results validate implementation and stated estimands; they
  do not establish biological epistasis or causality.

### Public real-genome subset validation

- A staged map → profile → associate → advanced-filter → report run completed
  on a public 12-genome SPN534 subset using real genome FASTA files and a
  documented alternative public reference.
- Nested 12/24/48-genome mapping and association runs provide engineering
  evidence. The alternative reference means these runs are not exact
  reproductions of the historical frozen SPN534 WGS mapping.

### User-deferred production-scale acceptance

- The complete 534-genome SPN534 reconstruction is deferred.
- The 224-genome × 25,000-locus production-scale run is deferred.
- The private Bacillus HPC validation is deferred.
- No production-readiness claim is made from microbenchmarks, controlled data,
  or the public SPN534 subsets.

## Installation

The indexed mapper is an optional dependency and Numba acceleration is a
separate optional dependency:

```bash
git clone --branch integration/unified-framework --single-branch \
  https://github.com/Naclist/Cinch.git
cd Cinch
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[mapping,performance]"
cinch --version
```

Install `.[test]` to run the complete suite. Install `.[release]` only when
building and checking distributions.

## Known limitations

- Historical `configure` and `uberBlast` are not redistributed; exact legacy
  mapper equivalence is unresolved.
- Indexed allele identifiers are deterministic sequence-hash identifiers and
  need not equal historical first-observation/uberBlast identifiers.
- Categorical conditional association currently supports one explicitly
  selected two-background contrast per run, not continuous covariates or a
  multi-background omnibus test.
- Conditional-exchangeability assumptions cannot be proven by software.
- Candidate-family q-values do not imply FDR control over hypotheses omitted by
  data-dependent pre-screening.
- Full 534-genome and production-scale acceptance remain deferred.
