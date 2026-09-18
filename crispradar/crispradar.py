"""
CRISPRadar: A Python library for scanning DNA sequences to identify candidate CRISPR guide RNAs.
"""

import argparse
import logging
import os
from Bio import SeqIO
import pandas as pd

try:
    from .scanner import scan_sequence
except ImportError:
    from scanner import scan_sequence


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scan DNA sequences to identify candidate CRISPR guide RNAs for a given PAM sequence."
    )

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        "-i", "--input_file",
        type=str,
        required=True,
        help="Path to the input FASTA file containing DNA sequence to scan."
    )

    # Optional arguments
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        "-l", "--guide_length",
        type=int,
        default=20,
        help="Length of the guide RNA (default: 20).",
    )
    optional.add_argument(
        "-p", "--pam",
        type=str,
        default="NGG",
        help="PAM sequence (default: NGG). Use N for any base.",
    )
    optional.add_argument(
        "-o", "--output_file",
        type=str,
        default="candidate_guides.csv",
        help="Path to the output CSV file (default: candidate_guides.csv).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # Check if input file exists
    if not os.path.isfile(args.input_file):
        logging.error(f"Input file '{args.input_file}' does not exist.")
        return 1

    # Scan sequences and collect candidate guides
    candidate_guides = []
    for record in SeqIO.parse(args.input_file, "fasta"):
        matches = scan_sequence(str(record.seq), args.guide_length, args.pam)
        candidate_guides.extend(
            {"record_id": record.id, **match.__dict__} for match in matches
        )

    columns = [
        "record_id",
        "guide",
        "pam",
        "guide_start",
        "guide_end",
        "pam_start",
        "pam_end",
        "strand",
    ]


    pd.DataFrame(candidate_guides, columns=columns).to_csv(
        args.output_file, index=False
    )
    logging.info("Wrote %d candidate guides to %s", len(candidate_guides), args.output_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
