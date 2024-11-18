import argparse
from Bio import SeqIO
import pandas as pd
import yaml

def load_config(config_file):
    """Load the configuration from a YAML file."""
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def main(args):
    # Load configuration
    config = load_config('config.yaml')
    num_probe_pairs = config.get('num_probe_pairs', 1)  # Default to 1 if not specified

    # Output file names based on config
    output_excel_combined = "output_after_cross_check_combined.xlsx"
    output_excel_lhs = "output_after_cross_check_lhs.xlsx"
    output_excel_rhs = "output_after_cross_check_rhs.xlsx"

    # Dictionary to store the LHS and RHS entries for each gene
    gene_entries = {}

    # Read the input FASTA file
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        gene_id = record.id.split("_")[0]
        entry_type = record.id.split("_")[3]
        index = record.id.split("_")[4] if len(record.id.split("_")) > 2 else ""

        # Initialize dictionary entry if it doesn't exist
        if gene_id not in gene_entries:
            gene_entries[gene_id] = {'LHS': [], 'RHS': []}

        # Append to the list of entries
        if entry_type in gene_entries[gene_id]:
            gene_entries[gene_id][entry_type].append((index, record))

    # Lists to store the data for Excel files
    combined_data = []
    lhs_data = []
    rhs_data = []

    # Collect the filtered entries for the Excel files
    for gene_id, entries in gene_entries.items():
        # Sort entries by their occurrence index
        lhs_records = sorted(entries.get("LHS", []), key=lambda r: r[0])
        rhs_records = sorted(entries.get("RHS", []), key=lambda r: r[0])

        # Create pairs up to num_probe_pairs
        num_pairs = min(len(lhs_records), len(rhs_records), num_probe_pairs)
        for i in range(num_pairs):
            lhs_record = lhs_records[i][1]
            rhs_record = rhs_records[i][1]
            combined_data.append([lhs_record.id, str(lhs_record.seq)])
            combined_data.append([rhs_record.id, str(rhs_record.seq)])
            lhs_data.append([lhs_record.id, str(lhs_record.seq)])
            rhs_data.append([rhs_record.id, str(rhs_record.seq)])

    # Create DataFrames and save to Excel files
    combined_df = pd.DataFrame(combined_data, columns=["Gene ID", "Sequence"])
    lhs_df = pd.DataFrame(lhs_data, columns=["Gene ID", "Sequence"])
    rhs_df = pd.DataFrame(rhs_data, columns=["Gene ID", "Sequence"])

    combined_df.to_excel(output_excel_combined, index=False)
    lhs_df.to_excel(output_excel_lhs, index=False)
    rhs_df.to_excel(output_excel_rhs, index=False)

    print(f"\nFiltered combined entries have been saved to {output_excel_combined}")
    print(f"Filtered LHS entries have been saved to {output_excel_lhs}")
    print(f"Filtered RHS entries have been saved to {output_excel_rhs}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Excel output based on filtered FASTA entries.")
    parser.add_argument("--input_fasta", required=True, help="Input FASTA file.")
    args = parser.parse_args()
    main(args)

