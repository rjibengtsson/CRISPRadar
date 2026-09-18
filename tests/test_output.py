from dataclasses import dataclass

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re
from dataclasses import dataclass


@dataclass(frozen=True)
class truthGuide:
    sequence: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class predictedGuide:
    sequence: str
    start: int
    end: int
    strand: str




def format_thruth_set(df):
    """
    Format the truth set dataframe to match the predicted dataframe format.
    """
    df = df.rename(columns={
        'guide_sequence': 'guide',
        'sstart': 'guide_start',
        'send': 'guide_end'
    })

    # Determine strand based on start and end positions
    df['strand'] = df.apply(lambda row: '-' if row['guide_start'] > row['guide_end'] else '+', axis=1)


    # Ensure start is less than end otherwise swap them
    df['guide_start'], df['guide_end'] = df[['guide_start', 'guide_end']].min(axis=1), df[['guide_start', 'guide_end']].max(axis=1)
    
    return df[['guide', 'guide_start', 'guide_end', 'strand']]



def format_predicted_set(df):
    """
    Format the predicted start value to match the truth set 1-index.
    """
    df['guide_start'] = df['guide_start'] + 1
    return df[['guide', 'guide_start', 'guide_end', 'strand']]




def populate_dataclass(df, dataclass_type):
    return [dataclass_type(row['guide'], row['guide_start'], row['guide_end'], row['strand']) for index, row in df.iterrows()]



def compare_guides(predicted_list, true_list):
    true_by_seq = {t.sequence: t for t in true_list}  # last one wins if duplicates
    rows = []

    for p in predicted_list:
        if p.sequence in true_by_seq and p.start == true_by_seq[p.sequence].start and p.end == true_by_seq[p.sequence].end:
            status = 'full_match'
        elif p.sequence not in true_by_seq:
            status = 'no_match'
        rows.append({'guide_sequence': p.sequence, 'strand': p.strand, 'match_status': status})

    return pd.DataFrame(rows)







# Truth set with location information 
df_truth = pd.read_csv("tests/data/truth_guides_blastn_results.txt", sep='\t', header=0)
df_truth_formatted = format_thruth_set(df_truth)
df_truth_dataclass = populate_dataclass(df_truth_formatted, truthGuide)


# Gather predicted set sequence list
df_predicted = pd.read_csv("tests/data/test.csv", header=0)
df_predicted_formated = format_predicted_set(df_predicted)
df_predicted_dataclass = populate_dataclass(df_predicted_formated, predictedGuide)


# Compare the two lists
matched_results = compare_guides(df_predicted_dataclass, df_truth_dataclass)

print(matched_results[matched_results['match_status'] == 'full_match'].shape[0], "full matches")
print(matched_results[matched_results['match_status'] == 'partial_match'].shape[0], "partial matches")
print(matched_results[matched_results['match_status'] == 'no_match'].shape[0], "no matches")
matched_results.to_csv("tests/data/CRISPRadar_performance_report.csv", index=False)


# print(matched_results[matched_results['match_status'] == 'full_match'].shape[0], "full matches")
# print(matched_results[matched_results['match_status'] == 'partial_match'].shape[0], "partial matches")
# print(matched_results[matched_results['match_status'] == 'no_match'].shape[0], "no matches")