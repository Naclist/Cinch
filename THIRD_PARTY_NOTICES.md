# Third-party notices and license scope

The MIT license in `LICENSE` applies to CINCH software and documentation owned
by Naclist, including source migrated from `Naclist/Cinch-dev2` into this
repository. It does not relicense third-party software, biological records,
publication data, or trademarks.

## Python dependencies

CINCH does not vendor its Python dependencies. `pip` obtains them separately,
and each remains governed by its own license.

| Dependency | Purpose | Upstream license |
|---|---|---|
| NumPy | arrays and numerical kernels | BSD-3-Clause and bundled component licenses |
| pandas | tabular data | BSD-3-Clause and bundled component licenses |
| PyArrow | Parquet/Arrow I/O | Apache-2.0 |
| SciPy | scientific calculations | BSD-3-Clause and bundled component licenses |
| scikit-learn | clustering/metrics | BSD-3-Clause |
| Matplotlib | report figures | PSF-based Matplotlib license and bundled component licenses |
| Biopython | sequence I/O | Biopython License Agreement |
| NetworkX | graph operations | BSD-3-Clause |
| PyYAML | YAML I/O | MIT |
| mappy/minimap2 | optional indexed mapping backend | MIT |
| Numba | optional compiled kernels | BSD-2-Clause and bundled component licenses |
| pytest | test runner | MIT |

Dependency license texts and notices are distributed by those projects, not
copied into the CINCH wheel. Release environments should retain the notices
installed with each dependency.

## Coinfinder and SPN534 evidence

CINCH contains no Coinfinder program source. Coinfinder itself is distributed
upstream under GPL-3.0. Compact CINCH development evidence under
`examples/spn534/` cites and is derived in part from the public
`fwhelan/coinfinder-manuscript` inputs and outputs associated with:

F. J. Whelan, M. Rusilowicz and J. O. McInerney (2020), "Coinfinder: detecting
significant associations and dissociations in pangenomes", *Microbial
Genomics* 6:e000338, DOI `10.1099/mgen.0.000338`.

The manuscript repository does not declare a repository-level software/data
license. Accordingly, `examples/spn534/` is excluded from the Python source and
wheel distributions and is not covered by CINCH's MIT grant. Its provenance is
recorded in `examples/spn534/PROVENANCE.md`. A future GitHub source release must
retain this notice and attribution; users who redistribute those evidence
artifacts are responsible for confirming the applicable upstream terms.

## NCBI records

The repository records versioned NCBI assembly accessions and SHA256 values but
does not redistribute the corresponding genome assemblies. NCBI records and
downloaded sequence files remain subject to NCBI's policies and the terms
attached to their submitters and source databases.

## External executables

Historical `configure` and `uberBlast` executables are not included. Their
redistribution and installation terms were not established and CINCH does not
grant rights to them.
