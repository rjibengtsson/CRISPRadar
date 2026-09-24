"""
CRISPRadar: A Python library for scanning DNA sequences to identify candidate CRISPR guide RNAs.
"""

import argparse
import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

try:
    from .scanner import scan_sequence
except ImportError:
    from scanner import scan_sequence


# Preset scanning parameters for supported Cas nucleases.
CAS_PRESETS = {
    "cas9": {"pam": "NGG", "pam_position": "3prime", "guide_length": 20},
    "cas12f": {"pam": "TTTR", "pam_position": "5prime", "guide_length": 20},
}
DEFAULT_CAS = "cas9"


def target_to_guide(target_seq: str) -> str:
    """
    Convert a DNA target (protospacer) sequence to its guide RNA spacer sequence.

    The target sequence is reported 5'->3' on the strand carrying the PAM, which
    is the same sequence and orientation as the guide RNA spacer, so transcribing
    it (T -> U) gives the guide RNA read 5'->3'.
    """
    return str(Seq(target_seq).transcribe())


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
        "-c", "--cas",
        type=str.lower,
        choices=sorted(CAS_PRESETS),
        default=DEFAULT_CAS,
        help="Cas nuclease preset: cas9 (NGG, 3prime, 20 nt) or cas12f "
             "(TTTR, 5prime, 20 nt). Explicit --pam, --pam_position and "
             f"--guide_length override the preset (default: {DEFAULT_CAS}).",
    )
    optional.add_argument(
        "-l", "--guide_length",
        type=int,
        default=None,
        help="Length of the guide RNA (default: from --cas preset).",
    )
    optional.add_argument(
        "-p", "--pam",
        type=str,
        default=None,
        help="PAM sequence. Use N for any base (default: from --cas preset).",
    )
    optional.add_argument(
        "--pam_position",
        type=str,
        choices=["5prime", "3prime"],
        default=None,
        help="Position of the PAM relative to the guide: 3prime (e.g. Cas9) "
             "or 5prime (e.g. Cas12a/Cas12f) (default: from --cas preset).",
    )
    optional.add_argument(
        "-o", "--output_file",
        type=str,
        default="candidate_guides.csv",
        help="Path to the output CSV file (default: candidate_guides.csv).",
    )
    args = parser.parse_args()

    # Fill any parameter not given explicitly from the selected Cas preset.
    for key, value in CAS_PRESETS[args.cas].items():
        if getattr(args, key) is None:
            setattr(args, key, value)
    return args


def main():
    args = parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # Check if input file exists
    if not os.path.isfile(args.input_file):
        logging.error(f"Input file '{args.input_file}' does not exist.")
        return 1

    logging.info(
        "Scanning with %s preset: PAM=%s, pam_position=%s, guide_length=%d",
        args.cas, args.pam, args.pam_position, args.guide_length,
    )

    # Scan sequences and collect candidate guides
    candidate_guides = []
    for record in SeqIO.parse(args.input_file, "fasta"):
        matches = scan_sequence(
            str(record.seq), args.guide_length, args.pam, pam_position=args.pam_position
        )
        candidate_guides.extend(
            {
                "record_id": record.id,
                **match.__dict__,
                "guide_seq": target_to_guide(match.target_seq),
            }
            for match in matches
        )

    columns = [
        "record_id",
        "target_seq",
        "guide_seq",
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
