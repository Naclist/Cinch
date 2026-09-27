# Performance benchmarks

## Production target

[COMPUTED | HIGH] `25,000 × 24,999 / 2 = 312,487,500` unordered locus pairs for 224 genomes on 20 CPU cores.

## Baseline evidence

| Component | Evidence | Result | Interpretation |
|---|---|---|---|
| frozen mapper threading | source inspection of `map_genomes`/`run_wgs` | serial; thread count recorded only | no claimed parallel speedup |
| frozen coordinate distance | source inspection | Python pair enumeration plus dense row for every pair | O(p²) rows and unsuitable at 25K |
| unified physical distance | controlled frozen comparison | sparse observed-pair source with exact reliable order/bp values | no dense missing-pair table; production memory unmeasured |
| frozen association | source inspection | nested Python pairs × four channels × DataFrame distance lookup | principal architectural risk |
| dev2 Numba kernels | source inspection and source tests | compiled triangular kernels exist in frozen scripts | candidate architecture; not unified benchmark evidence |
| unified weighted MI | 1,000,000 valid binary mass triplets, five repeats | best 0.083142542 s; 12,027,537 pairs/s | isolated vectorized kernel only; excludes state construction and I/O |
| unified indexed mapper | 12 synthetic genomes × 200 exact 300-nt loci | 1/2/4 workers: 0.534/0.560/0.468 s; 4-worker speedup 1.14× | small controlled workload; process/index overhead dominates; no biological claim |
| unified pair blocks | 224 samples × 80 loci; 3,160 pairs; four channels | 3.171 s and 2.99 MB traced peak versus frozen 3.212 s and 20.27 MB | exact row counts; 85.2% lower traced Python peak; I/O excluded |
| SHC permutation kernel | 1,000 samples; 20 SHCs; 199 permutations; one edge | 0.07095 s; 2,805 permutations/s | single run; excludes edge scheduling and I/O |
| sparse ARACNE | 500 nodes; 1,000 edges | 0.01542 s; 64,831 input edges/s | adjacency-triangle kernel; single run |
| public SPN map ladder | 12 real genomes × 50/100/200/400 alternative CDS | 2.56/2.66/2.88/3.83 s wall | four workers; includes process/import/index/output overhead |
| public SPN association ladder | 12 real genomes × 50/100/200/400 loci | 1.18/1.78/4.23/13.80 s wall | full PP/PT/TP/TT, before necessary-state prescreen |
| necessary-state prescreen | 12 real genomes × 200/400 loci | 4.23→2.41 s and 13.80→5.95 s | exact outputs; 1.76×/2.32× |
| optional compiled association | nested 12/24/48 real genomes × 400 alternative CDS; three repeats | Python medians 6.56/7.31/8.74 s; Numba medians 2.42/2.56/2.70 s | 2.71×/2.86×/3.24×; exact distance frames and non-floats; max float difference 2.22e-15 |
| public SPN48 E2E | 48 real genomes × 400 CDS, 99 permutations | map 7.31 s; associate 2.70 s median; filter 3.96 s; report 1.41 s | alternative reference; engineering/workflow evidence only |
| 25K synthetic mapper | 1 synthetic genome × 25,000 150-nt loci | 17.86 s wall; 302.4 MB RSS; 1,480 queries/s | engineering stress only |

## Benchmark ladder

| Level | Workload | Required measures | Status |
|---|---|---|---|
| 1 | synthetic correctness cases | wall, RSS, equality | COMPLETE: mapping controlled cases and exact pair A/B pass |
| 2 | small representative genomes/reference | mapping and pair throughput | COMPLETE: 12 public SPN genomes, 50–400 CDS |
| 3 | moderate subset | scaling and resume | PARTIAL: nested 12/24/48 sample ladder and resume; association core scaling absent |
| 4 | full reference, limited genomes | mapper/index memory | PARTIAL: synthetic 25K × 1; biological 25K reference absent |
| 5 | 224 genomes × 25K reference | mapping wall/CPU/RSS | BLOCKED: inputs absent |
| 6 | full pairwise workflow | pairs/s, wall/CPU/RSS/disk | BLOCKED: optional compiled scorer measured locally, but 25K all-pair resources remain unestimated |

## Compiled association A/B

[COMPUTED | HIGH] Nested 12/24/48-sample profiles from one 48-genome real mapping were run three times per engine over the same 79,800-pair universe. Python medians were 6.56/7.31/8.74 s and compiled medians were 2.42/2.56/2.70 s. Compiled throughput was 32,975/31,172/29,556 unordered pairs/s. Association remained effectively single-core, so this is a compiled-kernel gain rather than core scaling.

[COMPUTED | HIGH] At 48 samples, PP/PT/TP/TT row counts were 1,016/7,226/10,155/76,038 for both engines. At every sample size, complete physical-distance frames and all non-floating association fields were exact. The largest Python/Numba floating difference was `2.220446049250313e-15`, below the prespecified `1e-12` tolerance.

[COMPUTED | HIGH] At 48 samples, median maximum RSS decreased from 458,211,328 to 386,301,952 bytes (15.7%). At 12 samples it increased 2.4%, so the compiled path is not universally lower-memory. Evidence is stored in `benchmarks/results/spn534_public_subset_numba_20260927.json`.

[COMPUTED | HIGH] Profiling before compiled distance aggregation assigned 4.616 of 6.460 workflow seconds to sparse order/bp construction, 0.996 s to block scoring, and 0.406 s to Parquet writes. Compiling the sparse same-contig aggregation reduced the 48-sample total from a prior 5.24 s scorer-only median to 2.70 s while preserving the 66,621-row distance frame exactly.

## Optimization decision log ODL-001

Problem: frozen all-pair Python/Pandas construction is incompatible with 312,487,500 pairs.

Measured evidence: source complexity confirmed; no wall-time measurement yet.

Affected workflow: physical distance and PP/PT/TP/TT scoring.

Implemented intervention: deterministic blockwise triangular enumeration, necessary-state prescreening, precomputed feature states/cardinalities, an optional compiled four-channel contingency kernel, sparse observed distance lookup, and atomic block outputs.

Expected benefit: bounded memory and elimination of per-pair DataFrame allocation.

Scientific risk: altered eligibility masks, driver tie handling, ordering, or floating summation.

Implementation cost: high; the optional Numba dependency is isolated from the default install.

Validation plan: exhaustive small-case old/new comparison before moderate benchmark.

Stop condition: target block fits configured memory, exact pair/channel universe holds, and another stage dominates runtime. The local 400-locus gate passes, but the sample/core and 25K all-pair gates remain open.
