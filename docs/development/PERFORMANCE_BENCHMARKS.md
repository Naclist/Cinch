# Performance benchmarks

## Production target

[COMPUTED | HIGH] `25,000 × 24,999 / 2 = 312,487,500` unordered locus pairs for 224 genomes on 20 CPU cores.

## Baseline evidence

| Component | Evidence | Result | Interpretation |
|---|---|---|---|
| frozen mapper threading | source inspection of `map_genomes`/`run_wgs` | serial; thread count recorded only | no claimed parallel speedup |
| frozen coordinate distance | source inspection | Python pair enumeration plus dense row for every pair | O(p²) rows and unsuitable at 25K |
| frozen association | source inspection | nested Python pairs × four channels × DataFrame distance lookup | principal architectural risk |
| dev2 Numba kernels | source inspection and source tests | compiled triangular kernels exist in frozen scripts | candidate architecture; not unified benchmark evidence |
| unified weighted MI | 1,000,000 valid binary mass triplets, five repeats | best 0.083142542 s; 12,027,537 pairs/s | isolated vectorized kernel only; excludes state construction and I/O |
| unified indexed mapper | 12 synthetic genomes × 200 exact 300-nt loci | 1/2/4 workers: 0.534/0.560/0.468 s; 4-worker speedup 1.14× | small controlled workload; process/index overhead dominates; no biological claim |
| unified pair blocks | 224 samples × 80 loci; 3,160 pairs; four channels | 3.171 s and 2.99 MB traced peak versus frozen 3.212 s and 20.27 MB | exact row counts; 85.2% lower traced Python peak; I/O excluded |
| SHC permutation kernel | 1,000 samples; 20 SHCs; 199 permutations; one edge | 0.07095 s; 2,805 permutations/s | single run; excludes edge scheduling and I/O |
| sparse ARACNE | 500 nodes; 1,000 edges | 0.01542 s; 64,831 input edges/s | adjacency-triangle kernel; single run |

## Benchmark ladder

| Level | Workload | Required measures | Status |
|---|---|---|---|
| 1 | synthetic correctness cases | wall, RSS, equality | COMPLETE: mapping controlled cases and exact pair A/B pass |
| 2 | small representative genomes/reference | mapping and pair throughput | PARTIAL: synthetic mapper throughput measured; representative genomes absent |
| 3 | moderate subset | scaling and resume | NOT_STARTED |
| 4 | full reference, limited genomes | mapper/index memory | NOT_STARTED |
| 5 | 224 genomes × 25K reference | mapping wall/CPU/RSS | BLOCKED: inputs absent |
| 6 | full pairwise workflow | pairs/s, wall/CPU/RSS/disk | BLOCKED: streaming distance provider absent and resources unestimated |

## Optimization decision log ODL-001

Problem: frozen all-pair Python/Pandas construction is incompatible with 312,487,500 pairs.

Measured evidence: source complexity confirmed; no wall-time measurement yet.

Affected workflow: physical distance and PP/PT/TP/TT scoring.

Proposed intervention: deterministic blockwise triangular enumeration, precomputed feature states/marginals, compiled contingency kernels, sparse observed distance lookup, and atomic block outputs.

Expected benefit: bounded memory and elimination of per-pair DataFrame allocation.

Scientific risk: altered eligibility masks, driver tie handling, ordering, or floating summation.

Implementation cost: high.

Validation plan: exhaustive small-case old/new comparison before moderate benchmark.

Stop condition: target block fits configured memory, exact pair/channel universe holds, and another stage dominates runtime.
