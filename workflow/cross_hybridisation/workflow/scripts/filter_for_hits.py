import argparse
from Bio import SeqIO

def parse_blast_hits(blast_files, threshold=99.0):
    exclude_entries = set()
    for blast_file in blast_files:
        with open(blast_file, 'r') as f:
            for line in f:
                columns = line.split('\t')
                percent_identity = float(columns[2].strip())
                if percent_identity > threshold:
                    entry_name = columns[0].strip()
                    exclude_entries.add(entry_name)
    return exclude_entries

def filter_fasta(input_fasta, output_fasta, exclude_entries):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        write_seq = False
        for line in infile:
            if line.startswith('>'):
                entry_name = line[1:].strip()
                if entry_name not in exclude_entries:
                    outfile.write(line)
                    write_seq = True
                else:
                    write_seq = False
            else:
                if write_seq:
                    outfile.write(line)

def main(args):
    exclude_entries = parse_blast_hits(args.blast_outputs)
    filter_fasta(args.input_fasta, args.output_fasta, exclude_entries)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter FASTA entries based on BLAST hits.")
    parser.add_argument("--blast_outputs", required=True, nargs='+', help="BLAST output files.")
    parser.add_argument("--input_fasta", required=True, help="Input FASTA file.")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file.")
    args = parser.parse_args()
    main(args)

