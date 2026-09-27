# Conditional-association smoke example

This synthetic input checks the installed `profile` and
`conditional-associate` CLI paths. It is an installation smoke test, not a
power benchmark or biological validation dataset.

```bash
cinch profile examples/conditional_smoke/profile.tsv \
  -o conditional-profile --missing-policy absence

cinch conditional-associate \
  --profile conditional-profile/PROFILE_V2.npz \
  --metadata examples/conditional_smoke/metadata.tsv \
  --population-column population --background-column habitat \
  --reference-background soil --comparison-background hospital \
  --channels TT --pairs examples/conditional_smoke/pairs.tsv \
  --permutations 9 --seed 42 --min-eligible-samples 16 \
  --min-population-cell 4 --min-shared-populations 2 \
  --min-background-samples 8 -o conditional-results
```

The manifest must report `COMPLETE`, and the single TT differential hypothesis
must be tested. Nine permutations are deliberately insufficient for scientific
inference; they keep this installation check fast.
