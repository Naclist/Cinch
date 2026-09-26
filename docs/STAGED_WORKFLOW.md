# Staged unified workflow

## Available stages

```text
genome FASTA + reference CDS
        |
        v
cinch map --------------> PROFILE_V2.npz + coordinates + mapping manifest
        |                                      |
        |                                      v
legacy allele table --> cinch profile --> validated PROFILE_V2.npz
                                               |
                                               v
                                      cinch associate
                                               |
                                               v
                              distance table + restartable PP/PT/TP/TT blocks
                                               |
                                               v
                              cinch advanced-filter (PP binary only)
                                               |
                                               v
                                      cinch report
```

[KNOWN | HIGH] `cinch map` uses the indexed/resumable mapper. `cinch profile` converts an external allele table only when the user explicitly chooses whether non-calls mean absence or unresolved. `cinch associate` validates the profile, computes frozen order/bp distances, and writes deterministic restartable channel blocks.

## Commands

```bash
cinch map -r reference.cds.fasta -o mapped -t 4 genomes/*.fasta

cinch profile legacy.tsv -o converted --missing-policy unresolved

cinch associate \
  --profile mapped/profiles/PROFILE_V2.npz \
  --coordinates mapped/mapping/COORDINATES.tsv \
  -o associated --block-pairs 100000

cinch advanced-filter \
  --association associated/association \
  --profile mapped/profiles/PROFILE_V2.npz \
  --shc sample_shc.tsv --weights sample_weights.tsv \
  -o filtered --permutations 999

cinch report --filter-results filtered -o report
```

[KNOWN | HIGH] The frozen `cinch wgs` and `cinch filter` commands remain available and unchanged. The current `filter` consumes a complete frozen WGS result, not staged block output.

## Recovery boundary

[KNOWN | HIGH] Mapping caches each genome independently. Association hashes its state arrays, labels, distance table/schema, and configuration, then resumes only complete compatible channel blocks. Profile conversion and stage manifests use atomic replacement.

## Incomplete stage graph

[KNOWN | HIGH] `advanced-filter` applies the preserved binary SHC permutation, BH correction, and ARACNE only to PP edges. It does not misapply the binary test to PT/TP/TT categorical states. `report` is read-only with respect to statistical results.

[KNOWN | HIGH] Pair/channel scoring still enumerates the full hypothesis universe in Python, and the public SPN534 real-genome E2E validation remains pending. Therefore the staged workflow is not yet the production replacement for frozen `cinch wgs`/`filter`.
