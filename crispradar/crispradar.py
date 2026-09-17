"""
CRISPRadar: A Python library for scanning DNA sequences to identify candidate CRISPR guide RNAs.
"""

import argparse
import logging
import os, sys
from Bio import SeqIO
import pandas as pd
from scanner import CandidateGuide, reverse_complement, scan_sequence


