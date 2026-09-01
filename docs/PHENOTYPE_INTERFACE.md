# Phenotype interface

`--phenotype` accepts a delimited table whose first column contains unique WGS
sample IDs. Cinch v1 validates sample identity, freezes the exact bytes and
records SHA256. It does not test genotype-phenotype association, fit a GWAS,
select phenotype-specific edges or overlay phenotype on a network.

