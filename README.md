# CRISPRadar
Scans any input sequence to identify candidate CRISPR guide RNAs of a user-defined length adjacent to a specified PAM sequence.

## Python usage

```python
from crispradar import scan_sequence

matches = scan_sequence("AAAAAGGTTTTAGG", guide_length=4, pam="NRG")
```

`scan_sequence` returns every matching candidate guide on both strands. PAM sequences can use IUPAC degenerate bases such as `N`, `R`, and `Y`.
