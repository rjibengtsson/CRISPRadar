"""
This script is used to get the location of gRNA sequnece from the truth set (from PNAS paper and Scientific reports paper).
"""

import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path


def run_blastn_and_format(query_fasta, db_path, output_file, pident=85, cov_perc=80):
    blastn_cmd = [
        "blastn",
        "-task", "blastn-short",
        "-word_size", "7",
        # "-reward", "2",
        # "-penalty", "-3",
        # "-gapopen", "5",
        # "-gapextend", "2",
        # "-evalue", "0.05",
        "-query", f"{query_fasta}",
        "-db", f"{db_path}",
        "-max_target_seqs", "100",
        "-perc_identity", f"{pident}",
        "-qcov_hsp_perc", f"{cov_perc}",
        "-outfmt", "6",
        "-out", f"{output_file}"]
    print(f"Running BLASTN with command: {' '.join(blastn_cmd)}")
    subprocess.run(blastn_cmd, check=True)

    df = pd.read_csv(output_file, sep="\t", header=None,
                      names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    # df = df.loc[df.groupby("qseqid")["bitscore"].idxmax()]

    # map qseqid -> guide sequence from the query FASTA
    seq_map = {rec.id: str(rec.seq) for rec in SeqIO.parse(query_fasta, "fasta")}
    df["guide_sequence"] = df["qseqid"].map(seq_map)

    df.to_csv(output_file, sep="\t", index=False)
    return df




def generate_guide_fasta(df, output_dir):

    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    records = [
        SeqRecord(Seq(seq), id=str(seq_id), description="")
        for seq_id, seq in zip(df['No.'].astype(str) + '_' + df['Direction'], df["gRNA_sequence"])
    ]

    SeqIO.write(records, output_dir / "truth_guides.fasta", "fasta")
    return output_dir / "truth_guides.fasta"




def main():
    output_dir = Path("tests/data")

    # Load the truth set CSV file
    truth_set_csv = "tests/data/LTR_guide_sequences.csv"
    df = pd.read_csv(truth_set_csv, header=0)

    df = df.rename(columns={
        "Oligo ID": 'No.',
        "Sequence (5' to 3')": 'gRNA_sequence'
    })

    df = df[['Target name', 'Direction', 'No.', 'gRNA_sequence']]
    df["gRNA_sequence"] = df["gRNA_sequence"].str[4:]

    # for index, row in df.iterrows():
    #     if row['Direction'] == 'Sense':
    #         df.at[index, 'gRNA_sequence'] = str(Seq(row['gRNA_sequence'])[1:])
    #     elif row['Direction'] == 'Antisense':
    #         df.at[index, 'gRNA_sequence'] = str(Seq(row['gRNA_sequence'])[:-1])

    generate_guide_fasta(df, "tests/data")


    # Define the BLAST database path (make sure to create this database beforehand)
    db_path = "tests/data/AF324493.2_ltr"  # Replace with the path to your BLAST database

    query_fasta = output_dir / "truth_guides.fasta"
    blast_output_file = output_dir / "truth_guides_blastn_results.txt"

    # Run BLASTN to find locations of the guides in the target sequences
    run_blastn_and_format(query_fasta, db_path, blast_output_file)
    

if __name__ == "__main__":
    main()