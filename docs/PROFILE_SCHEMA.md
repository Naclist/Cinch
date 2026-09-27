# CINCH state-profile schema

## Contract

[KNOWN | HIGH] A CINCH profile has two aligned matrices with shape `samples × loci`.

| Matrix | Value | Meaning |
|---|---:|---|
| presence | `1` | locus detected as present |
| presence | `0` | locus confidently absent |
| presence | `-1` | presence unresolved/non-callable |
| types | `>=0` | callable nominal allele code |
| types | `-1` | allele non-callable |

[KNOWN | HIGH] A callable type is legal only when presence is `1`. Absence and unresolved states must have type `-1`. Allele codes are nominal identifiers, not ordered measurements.

## Legacy conversion

[KNOWN | HIGH] `from_legacy_alleles` requires `missing_policy="absence"` or `missing_policy="unresolved"`. It never guesses whether legacy zero, negative, `NA`, `NaN`, `None`, or blank values mean biological absence or technical non-callability.

[KNOWN | HIGH] Positive numeric values and non-missing textual allele labels are encoded independently within each locus. Code assignment is deterministic and the original value codebook is stored in metadata.

## Serialization

[KNOWN | HIGH] `save_profile` writes `samples`, `loci`, `presence`, `types`, and JSON metadata atomically. `load_profile` accepts both the V2 `types` key and frozen_v1 `type_state` key, then enforces the same invariants.

## Storage bound

[COMPUTED | HIGH] The two state matrices use five bytes per sample-locus cell: one byte for `int8` presence plus four bytes for `int32` type. At 224 samples × 25,000 loci, the raw matrices occupy 28,000,000 bytes (26.70 MiB), excluding labels, metadata, and temporary compression buffers.
