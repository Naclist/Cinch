# Release copyright and redistribution audit

Audit date: 2026-09-27
Target: Cinch v0.1.0 — Research Preview

## Decision

[COMPUTED | HIGH] Every commit in the integrated repository is authored by one
of the `Naclist` Git identities. The migrated Cinch-dev2 repository likewise
contains only `Naclist <776047333@qq.com>` commits. Both public source
repositories are owned by the GitHub account `Naclist`.

[KNOWN | HIGH] Before this release, `Naclist/Cinch` had no license and
`Naclist/Cinch-dev2` explicitly reserved all rights. The repository owner's
release instruction authorizes redistribution of the Naclist-owned integrated
software. The release therefore applies the MIT license prospectively to that
software and documentation while retaining source provenance.

[KNOWN | HIGH] This conclusion does not grant rights over third-party material.
The scope exceptions and dependency notices are recorded in
`THIRD_PARTY_NOTICES.md`.

## Integrated source audit

| Source | Integrated paths | Ownership/provenance result | Release treatment |
|---|---|---|---|
| `Naclist/Cinch` | `cinch/frozen_v1`, release documentation and assets | same repository owner/author identities | MIT |
| `Naclist/Cinch-dev2@ec1eaa5` | `cinch/statistics`, `cinch/population/shc.py`, `cinch/network/aracne.py`, `cinch/association/high_order.py` | same repository owner/author; migration manifest records each adaptation | MIT |
| preserved Cinch-dev2 snapshots | `cinch/legacy/dev2` | same owner/author; retained for scientific audit | MIT |
| CINCH Unified additions | remaining `cinch` modules, tests and development documentation | authored in this repository by Naclist | MIT |
| Coinfinder program | none | no source copied; comparator only | not distributed |
| `configure` / `uberBlast` | none | license/install contract unresolved | not distributed |

The file-level mapping is retained in `SOURCE_MIGRATION_MANIFEST.tsv`. No
copyright header or source-history attribution was removed.

## Direct dependency audit

The following packages are resolved from their own distributions and are not
vendored in CINCH. Their declared upstream licenses permit installation and use
alongside MIT-licensed CINCH; their own notice obligations remain separate.

| Dependency group | Packages | License result |
|---|---|---|
| core numerical/tabular | NumPy, pandas, SciPy, scikit-learn | BSD-family; NumPy/pandas also report bundled component licenses |
| storage | PyArrow | Apache-2.0 |
| plotting | Matplotlib | PSF-based Matplotlib license plus bundled component licenses |
| biology/graph/config | Biopython, NetworkX, PyYAML | Biopython agreement; BSD-3-Clause; MIT |
| optional mapping | mappy/minimap2 | MIT |
| optional performance | Numba | BSD-2-Clause plus bundled component licenses |
| test/release tooling | pytest, build, twine | permissive upstream licenses; development/release only |

The audit records direct dependencies, not every transitive dependency selected
by a future package resolver. Release environments must retain dependency
metadata/notices and should re-run the audit when dependency bounds change.

## External research data

[COMPUTED | HIGH] The Python sdist and wheel exclude `examples/spn534`. The
package artifacts therefore do not redistribute the public Coinfinder
manuscript inputs or the CINCH evidence derived from them.

[KNOWN | HIGH] The Git repository retains compact SPN534 evidence for audit
continuity. The upstream `fwhelan/coinfinder-manuscript` repository declares no
repository-level license. Those files are outside the CINCH MIT grant and carry
explicit provenance in `examples/spn534/PROVENANCE.md` and
`THIRD_PARTY_NOTICES.md`.

[KNOWN | HIGH] Versioned NCBI accessions and checksums are retained, while raw
assembly FASTA/GFF/GBFF files are excluded from source control and package
artifacts.

## Release gate

[KNOWN | HIGH] The MIT license is suitable for the Naclist-owned code because
no third-party program source is integrated and direct dependencies are linked
at installation rather than vendored. A release must still fail if artifact
inspection finds `examples/spn534`, raw genomes, `configure`, or `uberBlast` in
the sdist or wheel.
