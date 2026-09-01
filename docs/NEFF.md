# HC69 recurrence and Neff

For principal driver support counts `n_h` in HC69 background `h`:

```text
p_h  = n_h / sum_h(n_h)
Neff = 1 / sum_h(p_h^2)
```

Neff is the inverse Simpson concentration of background support. It equals one
when all support is in one background and rises as support spreads. It is not a
weighted-sample-size formula and is not an adjusted p-value.

The current frozen SPN534 rule first requires recurrence in at least three HC69
groups (`HC_supported >= 3`), then compares Neff with the order-matched background
mean. The later decision that future implementations may admit exactly two
populations (`HC_supported >= 2`) is a future specification only; the current
frozen result is not retroactively changed.

- `rsuA_2–glyS`: dominant TT driver restricted to one HC69; Neff=1; rejected.
- `gyrB–aguA`: support `7:1;16:3;27:3`; Neff=2.578947; retained.
