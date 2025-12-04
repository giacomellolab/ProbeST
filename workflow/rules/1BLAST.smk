# Rule to prepare BLAST database
rule prepare_blast_database:
    input:
        ref_genome="CDS_all.fa"
    output:
        expand("db/ref_genome-db.{suffix}", suffix=['ndb', 'nhr', 'nin', 'njs', 'nog', 'nos', 'not', 'nsq', 'ntf', 'nto'])
    log:
        "logs/prepare_blast_database.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        mkdir -p db
        makeblastdb -in {input.ref_genome} -dbtype nucl -parse_seqids -out db/ref_genome-db 2> {log}
        """

# Rule to run BLAST against reference
rule blast_against_ref:
    input:
        fasta="tmp_outputs/probes_hybparts_snakemake.fasta",
        db_files=rules.prepare_blast_database.output
    output:
        "tmp_outputs/probes_hybparts_off_targets.txt"
    params:
        db_dir="db",             # Directory where the database is located
        db="ref_genome-db",      # The name of the BLAST database
        evalue=30000,            # default 10, using 30,000 as value used in primer-blast
        word_size=7,             # Smaller word size for better sensitivity with short probes, default 11
        gapopen=3,               # Cost to open a gap
        penalty=-1,              # Penalty for mismatches to enforce specificity
        task="blastn-short"      # BLAST task optimized for short sequences, default blastn
    log:
        "logs/blast_against_ref.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        blastn -db {params.db_dir}/{params.db} -query {input.fasta} -out {output} \
               -outfmt "6 qseqid sseqid qstart qend sstart send pident mismatch" -evalue {params.evalue} -word_size {params.word_size} \
               -gapopen {params.gapopen} -penalty {params.penalty} -task {params.task} 2> {log}
        """

# Rule to filter probe pairs based on specificity (mismatches of each hit)
rule specificity_trim:
    input:
        "tmp_outputs/probes_hybparts_off_targets.txt"
    output:
        "tmp_outputs/trimmed_probes_hybparts_off_targets.txt"
    log:
        "logs/specificity_trim.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        python3 scripts/1BLAST/parse_blast_output.py {input} {output} 2> {log}
        """

# Rule to select final probe pairs
rule select_probes_pairs:
    input:
        primer3_output="tmp_outputs/primer3_output.txt",
        blast_output="tmp_outputs/trimmed_probes_hybparts_off_targets.txt"
    output:
        selected_probes="final_outputs/selected_probes.txt",
        probes_csv="final_outputs/probe_set.csv",
        probe_quantifications="final_outputs/probe_quantifications.txt",
        log_process = "final_outputs/process_log.txt" # internal log file of the different selection steps
    log:
        "logs/select_probes_pairs.log" # snakemake log file
    conda:
         "probes_env.yaml"
    shell:
        """
        mkdir -p logs
        mkdir -p final_outputs
        python3 scripts/1BLAST/select_probe_pairs.py {input.primer3_output} {input.blast_output} {output.probes_csv} {output.selected_probes} {output.probe_quantifications} {output.log_process} 2> {log}
        """
