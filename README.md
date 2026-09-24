# CRISPRadar

CRISPRadar is a lightweight Python toolkit for scanning DNA sequences for candidate CRISPR guide sites. It supports custom guide lengths, degenerate PAM sequences, and both forward and reverse-complement strand searching.

## Features

- Scan a DNA sequence for candidate guides matching a specified guide length and PAM
- Support degenerate PAM bases such as N, R, Y, S, W, K, M, B, D, H, and V
- Search both the forward strand and reverse complement
- Report 1-based inclusive coordinates for guide and PAM positions
- Export candidate guide hits as CSV from the command line
- Use directly from Python with a small, testable API

## Installation

From the project root:

```bash
python -m pip install biopython pandas numpy
```

## Command-line usage

The package includes a CLI for scanning a FASTA file and writing candidate guides to CSV.

```bash
python -m crispradar.crispradar -i tests/data/test.fasta -o candidate_guides.csv
```

Optional parameters:

```bash
python -m crispradar.crispradar \
  -i my_sequence.fasta \
  -c cas9 \
  -l 20 \
  -p NGG \
  --pam_position 3prime \
  -o candidate_guides.csv
```

Available presets:

- `cas9`: PAM `NGG`, PAM position `3prime`, guide length `20`
- `cas12f`: PAM `TTTR`, PAM position `5prime`, guide length `20`

The CLI writes a CSV with columns including:

- `record_id`
- `target_seq`
- `guide_seq`
- `pam`
- `guide_start`
- `guide_end`
- `pam_start`
- `pam_end`
- `strand`

## Python API

```python
from crispradar import scan_sequence

matches = scan_sequence(
    "AAAAAGGTTTTAGG",
    guide_length=4,
    pam="AGG",
    pam_position="3prime",
)

for match in matches:
    print(match)
```

This returns `CandidateGuide` objects with the matched guide, PAM, coordinates, and strand information.

```python
from crispradar import reverse_complement

print(reverse_complement("ATGRYN"))
# NRYCAT
```

## Example output

A typical match object looks like:

```python
CandidateGuide(
    target_seq='AAAA',
    pam='AGG',
    guide_start=1,
    guide_end=4,
    pam_start=5,
    pam_end=7,
    strand='+'
)
```

## Notes

- Sequences and PAMs must use supported nucleotide codes; unsupported bases raise a `ValueError`.
- PAM matching supports IUPAC ambiguity codes in the provided PAM sequence.
- Coordinates are reported in 1-based inclusive format, matching the convention used by the project tests.

## Repository layout

- `crispradar/scanner.py`: core scanning logic
- `crispradar/crispradar.py`: CLI and CSV export
- `tests/unit_test/test_scanner.py`: validation tests for guide scanning and coordinate behavior
