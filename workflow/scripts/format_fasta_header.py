'''
Sofia Rouot
Pyhton script to modify the FASTA headers of the S. Tm genome into the required format for the ProbeST pipeline.
'''
from Bio import SeqIO

input_file = "CDS_gene_targets.fa"
output_file = "CDS_gene_targets_formatted_headers.fasta"

with open(output_file, "w") as out_f:
    for record in SeqIO.parse(input_file, "fasta"):
        header = record.description

        # Extract fields from header
        gene = ""
        locus_tag = ""
        protein_id = ""

        if "[gene=" in header:
            gene = header.split("[gene=")[1].split("]")[0]
        if "[locus_tag=" in header:
            locus_tag = header.split("[locus_tag=")[1].split("]")[0]
        if "[protein_id=" in header:
            protein_id = header.split("[protein_id=")[1].split("]")[0]

        # Build new header
        new_header = f"gene_ID:{locus_tag} gene_name:{gene} transcript_ID:{protein_id}"
        record.description = new_header
        #record.id = new_header  # to replace the whole header line
        SeqIO.write(record, out_f, "fasta")