import csv
import pandas as pd 
import sys

def csv_to_fasta(csv_file, fasta_file):
    """
    Convert CSV -> FASTA with LHS (first flank_len bp) and RHS (last flank_len bp).
    Header format: >{id}_{gene}_{id}_LHS_{counter} and >{id}_{gene}_{id}_RHS_{counter}
    Counter is per (id, gene) pair and increments for each row.
    """
    counters = {}
    flank_len=25

    with open(csv_file, newline='') as inf, open(fasta_file, 'w') as outf:
        reader = csv.reader(inf)
        next(reader, None)
        
        for row in reader:
            if not row or len(row) < 3:
                continue  # skip bad lines

            seq_id = row[0].strip()
            seq = row[1].strip().upper()
            identifier = row[2].strip()

            # extract gene from identifier like: CBW16945.1|grxA|98e6b9c
            parts = identifier.split('|')
            gene = parts[1] if len(parts) > 1 else "unknown"

            # get LHS and RHS (if sequence shorter than flank_len this will just return what's available)
            lhs = seq[:flank_len]
            rhs = seq[-flank_len:] if len(seq) >= 1 else ""

            key = (seq_id, gene)
            count = counters.get(key, 0)

            # write LHS
            lhs_header = f">{seq_id}_{gene}_{seq_id}_LHS_{count}"
            outf.write(lhs_header + "\n" + lhs + "\n")

            # write RHS
            rhs_header = f">{seq_id}_{gene}_{seq_id}_RHS_{count}"
            outf.write(rhs_header + "\n" + rhs + "\n")

            # increment counter for this (id,gene)
            counters[key] = count + 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 prep_input_file.py <csv_file> <fasta_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    fasta_file = sys.argv[2]

    csv_to_fasta(csv_file, fasta_file)

