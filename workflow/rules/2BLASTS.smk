# Rule to prepare BLAST database for species 1 (pathogen)
rule prepare_blast_database_pathogen:
    input:
        ref_genome="CDS_all_pathogen.fa"
    output:
        expand("db_pathogen/ref_genome-db.{suffix}", suffix=['ndb', 'nhr', 'nin', 'njs', 'nog', 'nos', 'not', 'nsq', 'ntf', 'nto'])
    log:
        "logs/prepare_blast_database_pathogen.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        mkdir -p db_pathogen
        makeblastdb -in {input.ref_genome} -dbtype nucl -parse_seqids -out db_pathogen/ref_genome-db 2> {log}
        """

# Rule to run BLAST against reference for species 1 (pathogen)
rule blast_against_ref_pathogen:
    input:
        fasta="tmp_outputs/probes_hybparts_snakemake.fasta",
        db_files=rules.prepare_blast_database_pathogen.output
    output:
        "tmp_outputs/probes_hybparts_off_targets_pathogen.txt"
    params:
        db_dir="db_pathogen",    # Directory where the database is located
        db="ref_genome-db",      # The name of the BLAST database
        evalue=30000,            # default 10, using 30,000 as value used in primer-blast
        word_size=7,             # Smaller word size for better sensitivity with short probes, default 11
        gapopen=3,               # Cost to open a gap
        penalty=-1,              # Penalty for mismatches to enforce specificity
        task="blastn-short"      # BLAST task optimized for short sequences, default blastn
    log:
        "logs/blast_against_ref_pathogen.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        blastn -db {params.db_dir}/{params.db} -query {input.fasta} -out {output} \
               -outfmt "6 qseqid sseqid qstart qend sstart send pident mismatch" -evalue {params.evalue} -word_size {params.word_size} \
               -gapopen {params.gapopen} -penalty {params.penalty} -task {params.task} 2> {log}
        """

# Rule to prepare BLAST database for species 2 (host)
rule prepare_blast_database_host:
    input:
        ref_host="CDS_all_host.fa"
    output:
        expand("db_host/ref_host-db.{suffix}", suffix=['ndb', 'nhr', 'nin', 'njs', 'nog', 'nos', 'not', 'nsq', 'ntf', 'nto'])
    log:
        "logs/prepare_blast_database_host.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        mkdir -p db_host
        makeblastdb -in {input.ref_host} -dbtype nucl -parse_seqids -out db_host/ref_host-db 2> {log}
        """


# Rule to run BLAST against reference for species 2 (host)
rule blast_against_ref_host:
    input:
        fasta="tmp_outputs/probes_hybparts_snakemake.fasta",
        db_files=rules.prepare_blast_database_host.output
    output:
        "tmp_outputs/probes_hybparts_off_targets_host.txt"
    params:
        db_dir="db_host",        # Directory where the database is located
        db="ref_host-db",      # The name of the BLAST database
        evalue=30000,            # default 10, using 30,000 as value used in primer-blast
        word_size=7,             # Smaller word size for better sensitivity with short probes, default 11
        gapopen=3,               # Cost to open a gap
        penalty=-1,              # Penalty for mismatches to enforce specificity
        task="blastn-short"      # BLAST task optimized for short sequences, default blastn
    log:
        "logs/blast_against_ref_host.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        blastn -db {params.db_dir}/{params.db} -query {input.fasta} -out {output} \
               -outfmt "6 qseqid sseqid qstart qend sstart send pident mismatch" -evalue {params.evalue} -word_size {params.word_size} \
               -gapopen {params.gapopen} -penalty {params.penalty} -task {params.task} 2> {log}
        """

# Rule to filter probe pairs based on specificity for pathogen CDS and host CDS
rule specificity_trim_pat_host:
    input:
        pat_blast="tmp_outputs/probes_hybparts_off_targets_pathogen.txt",
        host_blast="tmp_outputs/probes_hybparts_off_targets_host.txt"
    output:
        trimmed_pat="tmp_outputs/trimmed_probes_hybparts_off_targets_pat.txt",
        trimmed_host="tmp_outputs/trimmed_probes_hybparts_off_targets_host.txt"
    log:
        "logs/specificity_trim_pat_host.log"
    conda:
         "probes_env.yaml"
    shell:
        """
        python3 scripts/2BLASTS/target_specificity_trim_pat_host.py {input.pat_blast} {input.host_blast} {output.trimmed_pat} {output.trimmed_host} > {log} 2>&1
        """

# Rule to select final probe pairs
rule select_probes_pairs_two_blast_rounds:
    input:
        primer3_output="tmp_outputs/primer3_output.txt",
        pat_blast_output="tmp_outputs/trimmed_probes_hybparts_off_targets_pat.txt",
        host_blast_output="tmp_outputs/trimmed_probes_hybparts_off_targets_host.txt"
    output:
        selected_probes="final_outputs/selected_probes.txt",
        probes_csv="final_outputs/probe_set.csv",
        probe_quantifications="final_outputs/probe_quantifications.txt",
        log_process = "final_outputs/process_log.txt" # internal log file of the different selection steps
    log:
        "logs/select_probes_pairs_two_blast_rounds.log" # snakemake log file
    conda:
         "probes_env.yaml"
    shell:
        """
        mkdir -p logs
        python3 scripts/2BLASTS/select_probe_pairs_two_blast_rounds_v2.py {input.primer3_output} {input.pat_blast_output} {input.host_blast_output} {output.probes_csv} {output.selected_probes} {output.probe_quantifications} {output.log_process} 2> {log}
        """
