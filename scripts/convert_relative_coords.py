#!/usr/bin/env python3
"""
Convert blastn -outfmt 6 coordinates to positions relative to a reference
point (e.g. R start of the LTR), using the convention:
    reference position -> +1
    one base upstream   -> -1
    (no position 0)

Usage:
    python convert_to_relative_coords.py blast_results.tsv 454 output.tsv
"""

import sys
import pandas as pd
import numpy as np

# Default blastn -outfmt 6 column names
DEFAULT_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "guide_sequence"
]

# Columns to convert (edit if your qstart/sstart are on a different
# coordinate system, e.g. only sstart/send are on the genome you care about)
COLS_TO_CONVERT = ["sstart", "send"]


def to_relative_array(pos: np.ndarray, ref: int) -> np.ndarray:
    """Vectorised conversion of absolute positions (numpy array) to relative positions."""
    if pos >= ref:
        return pos - (ref - 1)
    else:
        return pos - ref



def main(blast_file: str, ref_pos: int, out_file: str) -> None:
    # Assumes standard tab-separated -outfmt 6 with no header.
    # If your file has a header or custom columns, adjust `names=` / `header=`.
    df = pd.read_csv(blast_file, sep="\t", header=0)



    for col in COLS_TO_CONVERT:
        df[col] = df[col].astype(int)  # Ensure positions are integers
        df[col] = to_relative_array(df[col].to_numpy(), ref_pos)

    df.to_csv(out_file, sep="\t", index=False)
    print(f"Converted {len(df)} rows. Reference position {ref_pos} -> +1. Saved to {out_file}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <blast_outfmt6.tsv> <reference_absolute_position> <output.tsv>")
        sys.exit(1)

    blast_input = sys.argv[1]
    reference_position = int(sys.argv[2])
    output_path = sys.argv[3]

    main(blast_input, reference_position, output_path)