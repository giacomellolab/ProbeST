"""
Snakemake pipeline

This file is used to generate a txt file to be used for Primer3 in the command line / next rule in Snakefile, for designing probes for 
(a selection of) genes, according to the probe guidelines as provided by 10X genomics.

Input:
- the CDS fasta file which should be named CDS_gene_targets.fa. Contains only headers and sequences of genes of interest. MUST DO: user must manually modify headers for the following format:
        - position 0 --> gene ID
        - position 1 --> gene name
        - position 2 --> transcript name

Note: The headers of the fasta file should at least contain "gene_id"|"gene_name|"transcript_id". It can contain additional information for the user's own use. 

The sequence ID in the fasta and the primer3 input file is a "_" separated ID of gene ID, gene name and transcript ID.


Output:
- primer3_input_snakemake.txt. A txt file containing primer3 instructions.

"""

# Load packages  
import sys
import pandas as pd
from Bio import SeqIO 

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise ValueError("Usage: python primer3_input_design_Snakemake_pipeline.py CDS_gene_targets.fa")

    genes_of_interest = sys.argv[1]


    try:
        seq_df = pd.DataFrame(columns = ["Sequence", "gene_id", "gene_name", "transcript_id"])
        for record in SeqIO.parse(genes_of_interest, "fasta"):
            Sequence = str(record.seq)
            gene_id = record.description.split(" ")[0]
            gene_id = gene_id.split(":")[1]
            gene_name = record.description.split(" ")[1]
            gene_name = gene_name.split(":")[1]
            transcript_id = record.description.split(" ")[2]
            transcript_id = transcript_id.split(":")[1]
            df_temp = pd.DataFrame({"Sequence": Sequence, "gene_id": gene_id, "gene_name": gene_name, "transcript_id": transcript_id}, index = [0])
            seq_df = pd.concat([df_temp, seq_df.loc[:]]).reset_index(drop=True)

        # Write primer3 input file
        output_file = "primer3_input_snakemake.txt"
        with open(output_file, 'w') as f:
            for index, row in seq_df.iterrows():
                f.write('SEQUENCE_ID=' + row['gene_id'] + "_" + row['gene_name'] + "_" + row['transcript_id'] + '\n' + 'SEQUENCE_TEMPLATE=' + row['Sequence'] + '\n')
                f.write('PRIMER_PICK_LEFT_PRIMER=0' + '\n') #indicate creation of left side probes (identical to CDS)
                f.write('PRIMER_PICK_INTERNAL_OLIGO=0' + '\n') #indicate creation of interal hybridising probes 
                f.write('PRIMER_PICK_RIGHT_PRIMER=1' + '\n') #indicate creation of right side probes (reverse complement to CDS) 
                f.write('PRIMER_MIN_SIZE=25' + '\n') #probe size
                f.write('PRIMER_OPT_SIZE=25' + '\n') #probe size
                f.write('PRIMER_MAX_SIZE=25' + '\n') #probe size
                f.write('PRIMER_MAX_NS_ACCEPTED=1' + '\n') #max number of unknown nucleotides 
                f.write('PRIMER_RIGHT_MIN_GC=44' + '\n') #minimum GC percentage
                f.write('PRIMER_RIGHT_MAX_GC=72' + '\n') #maximum GC percentage
                f.write('PRIMER_RIGHT_MAX_POLY_X=5' + '\n') #maximum length of mononucleotide repeats
                f.write('PRIMER_NUM_RETURN=100' + '\n') #number of primers to return
                f.write('PRIMER_EXPLAIN_FLAG=1' + '\n') #get statistics of primer
                f.write('=' + '\n') #indicate end of parameters

    except Exception as e:
        print("An error occurred:", str(e))
