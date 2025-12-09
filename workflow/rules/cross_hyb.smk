configfile: "config/cross-hyb.yaml"

# Add your directory of blast_db_dir in the configfile.

# Access configuration variables
BLAST_DB_DIR = config["blast_db_dir"]
BLAST_DBS = config["blast_dbs"]

# Rule to prepare probes list from final main output
rule prepare_probeslist:
    input:
        SELECTED_PROBES
    output:
        "tmp_outputs/cross_hyb_outputs/probes.fasta"
    log:
        "logs/cross_hyb_logs/prepare_probeslist_cross_hyb.log"
    conda:
         "probes_env.yaml"
    benchmark:
        "benchmarks/rule_prepare_probeslist.txt"
    shell:
        """
        mkdir -p logs/cross_hyb_logs
	mkdir -p tmp_outputs/cross_hyb_outputs
        python3 scripts/cross-hyb/extract_probes_from_part1.py {input} {output} 2> {log}
        """

# Rule to merge LHS and RHS for each gene
rule merge_LHS_RHS_for_each_gene:
    input:
        "tmp_outputs/cross_hyb_outputs/probes.fasta"
    output:
        "tmp_outputs/cross_hyb_outputs/merged_sequences.txt"
    log:
        "logs/cross_hyb_logs/merge_LHS_RHS_for_each_gene.log"
    conda:
         "probes_env.yaml"
    benchmark:
        "benchmarks/rule_merge_LHS_RHS_for_each_gene.txt"
    shell:
        """
        python3 scripts/cross-hyb/script_merge_probes.py {input} {output} 2> {log}
        """

# Rule to trim from A-tail for cross-hybridization check
rule trim_from_Atail_for_cross_hybridization_check:
    input:
        "tmp_outputs/cross_hyb_outputs/merged_sequences.txt"
    output:
        "tmp_outputs/cross_hyb_outputs/cleaned_sequences.fasta"
    log:
        "logs/cross_hyb_logs/trim_from_Atail_for_cross_hybridization_check.log"
    conda:
         "probes_env.yaml"
    benchmark:
        "benchmarks/rule_trim_from_Atail_for_cross_hybridization_check.txt"
    shell:
        """
        python3 scripts/cross-hyb/clean_from_Atail.py {input} {output} 2> {log}
        """

# Rule to trim from prefix for cross-hybridization check
rule trim_from_prefix_for_cross_hybridization_check:
    input:
        "tmp_outputs/cross_hyb_outputs/cleaned_sequences.fasta"
    output:
        "tmp_outputs/cross_hyb_outputs/cleaned_sequences_from_prefix.fasta"
    log:
        "logs/cross_hyb_logs/trim_from_prefix_for_cross_hybridization_check.log"
    conda:
         "probes_env.yaml"
    benchmark:
        "benchmarks/rule_trim_from_prefix_for_cross_hybridization_check.txt"
    shell:
        """
        python3 scripts/cross-hyb/clean_prefix.py {input} {output} 2> {log}
        """

# Rule to download the BLAST database if not already downloaded
#rule download_db_blast_for_cross_hybridization_check:
#    output:
#        "download_complete.txt"
#    run:
#        import os
#        if not os.path.exists(output[0]):
#            shell(f"""
#                chmod +x download_db.sh
#                ./download_db.sh {BLAST_DB_DIR}
#                touch {output[0]}
#            """)
#        else:
#             print("Database already downloaded, skipping download step.")

rule blast_for_cross_hybridization_check:
    input:
        fasta="tmp_outputs/cross_hyb_outputs/cleaned_sequences_from_prefix.fasta"
    output:
        expand("tmp_outputs/blast_cross_hyb_results/results_{db}_nt.txt", db=BLAST_DBS)
    params:
        db_dir=BLAST_DB_DIR,
        dbs=BLAST_DBS,
        evalue=1e-10,  # Increase e-value for lower sensitivity
        word_size=15,  # Increase word size for lower sensitivity
        reward=1,      # Default is 1, lower reward might decrease sensitivity
        penalty=-1,    # Default is -2, increase penalty might decrease sensitivity
        task="blastn-short"  # Task optimized for short sequences
    log:
        "logs/cross_hyb_logs/blast_for_cross_hybridization_check.log"
    conda:
         "probes_env.yaml"
    benchmark:
        "benchmarks/rule_blast_for_cross_hybridization_check.txt"
    shell:
        """
        mkdir -p tmp_outputs/blast_cross_hyb_results
        for db in {params.dbs}; do
            blastn -db {params.db_dir}/$db -query {input.fasta} -out tmp_outputs/blast_cross_hyb_results/results_${{db}}_nt.txt -outfmt 6 -evalue {params.evalue} -word_size {params.word_size} -reward {params.reward} -penalty {params.penalty} -task {params.task} 2>> {log}
        done
        """

# Rule to filter hits from cross-hybridization check
# If no hits, the input and output will look the same
rule filter_hits:
    input:
        blast_outputs=expand("tmp_outputs/blast_cross_hyb_results/results_{db}_nt.txt", db=BLAST_DBS),
        fasta="tmp_outputs/cross_hyb_outputs/probes.fasta"
    output:
        "tmp_outputs/cross_hyb_outputs/filtered_probes.fasta"
    log:
        "logs/cross_hyb_logs/filter_hits.log"
    conda:
         "probes_env.yaml"
    benchmark:
        "benchmarks/rule_filter_hits.txt"
    shell:
        """
        python scripts/cross-hyb/filter_for_hits.py --blast_outputs {input.blast_outputs} --input_fasta {input.fasta} --output_fasta {output} 2> {log}
        """

rule generate_output_after_cross_check:
    input:
        primer3_output="tmp_outputs/primer3_output.txt",
	fasta="tmp_outputs/cross_hyb_outputs/filtered_probes.fasta"
    output:
        selected_probes="final_outputs/cross_hyb_filter/selected_probes.txt",
        probes_csv="final_outputs/cross_hyb_filter/probe_set.csv",
        probe_quantifications="final_outputs/cross_hyb_filter/probe_quantifications.txt",
        log_process = "final_outputs/cross_hyb_filter/process_log.txt" # internal log file of the different selection steps
    log:
        "logs/cross_hyb_logs/generate_output_after_cross_check.log"
    conda:
         "probes_env.yaml"
    benchmark:
        "benchmarks/rule_generate_output_after_cross_check.txt"
    shell:
        """
	mkdir -p final_outputs/cross_hyb_filter
        python scripts/cross-hyb/generate_output_with_hits.py {input.primer3_output} {input.fasta} {output.probes_csv} {output.selected_probes} {output.probe_quantifications} {output.log_process} 2> {log}
        """
