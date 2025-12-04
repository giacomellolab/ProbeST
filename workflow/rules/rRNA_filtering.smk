# Step to filter for rRNA from SILVA database

# Rule to prep input file probe_set.csv into fasta file
rule prep_input_file:
    input:
        csv_file="final_outputs/probe_set.csv"
    output:
        fasta_file="tmp_outputs/rRNA/probe_set.fasta"
    log:
        "logs/rRNA_filtering/prep_input_file.log"
    conda:
         "probes_env.yaml"
    shell:
        """
	mkdir -p logs/rRNA_filtering
	mkdir -p final_outputs/rRNA_filtering
	mkdir -p tmp_outputs/rRNA
        python3 scripts/rRNA/prep_input_file.py {input.csv_file} {output.fasta_file} > {log} 2>&1
        """

# Rule to prepare BLAST database for rRNA
rule prepare_blast_database_rRNA:
    input:
        ref_genome="SILVA_138.2_SSURef_NR99_tax_silva.fasta"
    output:
        expand("db_rRNA/rRNA-db.{suffix}", suffix=['ndb', 'nhr', 'nin', 'njs', 'nog', 'nos', 'not', 'nsq', 'ntf', 'nto'])
    log:
        "logs/rRNA_filtering/prepare_blast_database_rRNA.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        mkdir -p db_rRNA
        makeblastdb -in {input.ref_genome} -dbtype nucl -parse_seqids -out db_rRNA/rRNA-db 2> {log}
        """

# Rule to run BLAST against rRNA db
rule blast_against_rRNA_db:
    input:
        fasta="tmp_outputs/rRNA/probe_set.fasta",
        db_files=rules.prepare_blast_database_rRNA.output
    output:
        "tmp_outputs/rRNA/probes_off_targets_rRNA.txt"
    params:
        db_dir="db_rRNA",    # Directory where the database is located
        db="rRNA-db",      # The name of the BLAST database
        evalue=30000,            # 30,000, testing 1e-5
        word_size=7,             # Smaller word size for better sensitivity with short probes, default 11
        gapopen=3,               # Cost to open a gap
        penalty=-1,              # Penalty for mismatches to enforce specificity
        task="blastn-short"      # BLAST task optimized for short sequences, default blastn
    log:
        "logs/rRNA_filtering/blast_against_rRNA_db.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        blastn -db {params.db_dir}/{params.db} -query {input.fasta} -out {output} \
               -outfmt "6 qseqid sseqid qstart qend sstart send pident mismatch" -evalue {params.evalue} -word_size {params.word_size} \
               -gapopen {params.gapopen} -penalty {params.penalty} -task {params.task} 2> {log}
        """


# Rule to filter probe pairs based on specificity
rule specificity_trim_rRNA:
    input:
        blast_output="tmp_outputs/rRNA/probes_off_targets_rRNA.txt"
    output:
        trimmed_blast_output="tmp_outputs/rRNA/trimmed_probes_off_targets_rRNA.txt"
    log:
        "logs/rRNA_filtering/specificity_trim_rRNA.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        python3 scripts/rRNA/target_specificity_trim_rRNA.py {input.blast_output} {output.trimmed_blast_output} > {log} 2>&1
        """

# Rule to select final probe pairs after rRNA off-target filtering
rule select_probes_pairs_rRNA:
    input:
        primer3_output="tmp_outputs/primer3_output.txt",
        trimmed_blast_output="tmp_outputs/rRNA/trimmed_probes_off_targets_rRNA.txt"
    output:
        selected_probes="final_outputs/rRNA_filtering/selected_probes_rRNA.txt",
        probes_csv="final_outputs/rRNA_filtering/probe_set_rRNA.csv",
        probe_quantifications="final_outputs/rRNA_filtering/probe_quantifications_rRNA.txt",
        log_process = "final_outputs/rRNA_filtering/process_log_rRNA.txt" # internal log file of the different selection steps
    log:
        "logs/rRNA_filtering/select_probes_pairs_rRNA.log" # snakemake log file
    conda:
         "probes_env.yaml"
    shell:
        """
        mkdir -p logs
        python3 scripts/rRNA/select_probe_pairs_rRNA.py {input.primer3_output} {input.trimmed_blast_output} {output.probes_csv} {output.selected_probes} {output.probe_quantifications} {output.log_process} 2> {log}
        """

