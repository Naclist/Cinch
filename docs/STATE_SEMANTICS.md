# State semantics

`P` is `-1` noncallable, `0` confidently absent or `1` present. `T` is a
non-negative nominal, locus-local complete-CDS type only when `P=1`; otherwise
it is `-1`. Present ambiguous multicopy or incomplete type calls retain `P=1,
T=-1`. Missing and absence are never merged.

PP uses samples callable at both presence loci. PT requires B present and typed;
TP requires A present and typed; TT requires both loci present and typed. TT is
not restricted to core loci and uses pairwise, not global, completeness.

