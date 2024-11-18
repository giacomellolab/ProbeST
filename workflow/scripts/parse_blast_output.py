"""
Snakemake pipeline. Author: Sofia Rouot.

Since BLAST performs local alignment, the output includes matches which are not the full query length
    - The matches which are less than 20 nucleotides long are deleted. 
    - For the matches which are between 20 and 24 nucleotides, we correct the number of mismatches to represent full query alignment. 
    - The new column "n_corrected_mismatches" represents the amount of mismatches for full query alignment. 

The script also removes duplicate hits of a probe against multiple transcript variants of the same gene. This would occur when using a transcriptome database, for instance. The hit with the highest percentage identity score to the probe is kept.

Hits with over 6 mistmatches are removed, as these do not represent potential off-target hybridisation events. Over 6 mismatches is considered as enough to prevent off-target hybridisation.

"""

import pandas as pd

# Read blast output as pandas dataframe 
df = pd.read_csv("probes_hybparts_off_targets.txt", 
                 sep = "\t",
                 names = ["probe_id", "transcript_id", "qstart", "qend", "sstart", "send", "percentage_ident", "n_mismatches"]) 

# Add a column for the corrected mismatches, 20 is an arbitrary number
df['n_corrected_mismatches'] = "20"

# Convert transcript_id column to string and split to extract gene_id
df['gene_id'] = df['transcript_id'].astype(str).str.split('.').str[0]

# Remove duplicate combinations of probe_id and gene_id.
df = df.sort_values(["probe_id", "gene_id", "n_mismatches"]) #sort first on probe_id, then gene_id and then percentage identity
# Drop duplicates based on probe_id and gene_id, keeping the first occurrence
df = df.drop_duplicates(subset=['probe_id', 'gene_id'], keep='first').reset_index(drop=True)


# Trim blast output to only contain full length alignments
# Calculate the alignment length column
df['alignment_length'] = df['qend'] - df['qstart'] + 1

# Filter out rows where the alignment length is less than 20 nucleotides
df = df[df['alignment_length'] >= 20]

# Calculate the corrected mismatches for alignments between 20 and 24 nucleotides
mask_20_24 = (df['alignment_length'] >= 20) & (df['alignment_length'] <= 24)
df.loc[mask_20_24, 'n_corrected_mismatches'] = df['n_mismatches'] + (25 - df['alignment_length'])

# For alignments of exactly 25 nucleotides, the corrected mismatches are the same as 'n_mismatches'
df.loc[df['alignment_length'] == 25, 'n_corrected_mismatches'] = df['n_mismatches']

# Check for alignments longer than 25 nucleotides and raise an error or handle it
if (df['alignment_length'] > 25).any():
    print("Error: Probe is longer than 25 nucleotides")


# With the corrected mismatch score, delete hits which have more than 5 mismatches. 
df = df[df['n_corrected_mismatches'] < 6]

# Sort dataframe
df = df.sort_values(["probe_id", "n_corrected_mismatches"]) #sort form high to low corrected percentage identity 

# Export datafarme to csv 
df.to_csv("trimmed_probes_hybparts_off_targets.txt", sep=" ", index=False, header=True)

