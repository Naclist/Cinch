# Coinfinder concordance benchmark

The comparator audit contains 985 high-MI candidate pairs: 582 match published
Coinfinder association edges, 400 match dissociation edges, 982/985 (99.70%) are
therefore direct published edges, and 963/985 (97.77%) occur in the same published
component. No candidate lies wholly outside published components.

This establishes technical concordance with a published dependency result, not
sensitivity/specificity against biological truth. SPN534 was inspected during
development and is not an unseen test.

The 30-pair V-ATPase control shows why detection and linkage interpretation must
remain separate: strong dependency should be detected, then labeled local-linked.
Exact tables are under `examples/spn534/diagnostics`.
