import argparse
#from Bio import SeqIO
import pandas as pd


def parse_blast_hits(blast_files):
    """
    Parse BLAST outputs and return probe_ids that have
    ≥1 significant hit (to be excluded).
    """

    significant_hits = []

    for blast_file in blast_files:
        df = pd.read_csv(
            blast_file,
            sep="\t",
            names=[
                "probe_id", "transcript_id",
                "qstart", "qend",
                "sstart", "send",
                "percentage_ident", "n_mismatches"
            ]
        )

        # Extract gene_id
        df["gene_id"] = df["transcript_id"].astype(str).str.split(".").str[0]

        # Alignment length
        df["alignment_length"] = df["qend"] - df["qstart"] + 1

        # Keep full-length alignments
        df = df[df["alignment_length"] >= 20]

        # Correct mismatches
        df["n_corrected_mismatches"] = df["n_mismatches"]
        mask_20_24 = df["alignment_length"].between(20, 24)
        df.loc[mask_20_24, "n_corrected_mismatches"] += (
            25 - df.loc[mask_20_24, "alignment_length"]
        )

        # Reject weak hits
        df = df[df["n_corrected_mismatches"] < 6]

        # Keep best hit per probe × gene
        df = (
            df.sort_values(["probe_id", "gene_id", "n_corrected_mismatches"])
              .drop_duplicates(subset=["probe_id", "gene_id"], keep="first")
        )

        significant_hits.append(df)

    # Combine all BLAST outputs
    hits = pd.concat(significant_hits, ignore_index=True)

    # Any probe appearing at least once is non-specific
    exclude_entries = set(hits["probe_id"].unique())

    return exclude_entries


def filter_fasta(input_fasta, output_fasta, exclude_entries):
    """
    Remove FASTA entries whose headers match excluded probe IDs.
    """
    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        write_seq = False

        for line in infile:
            if line.startswith(">"):
                entry_name = line[1:].strip()
                write_seq = entry_name not in exclude_entries
                if write_seq:
                    outfile.write(line)
            else:
                if write_seq:
                    outfile.write(line)


def main(args):
    exclude_entries = parse_blast_hits(args.blast_outputs)
    filter_fasta(args.input_fasta, args.output_fasta, exclude_entries)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter FASTA probes using corrected-mismatch BLAST specificity"
    )
    parser.add_argument(
        "--blast_outputs",
        required=True,
        nargs="+",
        help="BLAST output files (tab-separated)"
    )
    parser.add_argument(
        "--input_fasta",
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "--output_fasta",
        required=True,
        help="Filtered FASTA file"
    )

    args = parser.parse_args()
    main(args)
