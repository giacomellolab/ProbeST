# ProbeST: a custom probe design pipeline for sequencing-based spatial transcriptomics 

Sofia Rouot, Ireen van Dolderen, Patrick Rosendahl Andreassen, Solène Frapard,  Sybil A. Herrera-Foessel, Hailey Sounart, Sami Saarenpää, Julia Vorholt, Stefania Giacomello

ProbeST allows to custom design binding probes for any gene of any species to be used in probe-based spatial transcriptomics protocols. To date, the Visium platform provides probes for the human and mouse transcriptomes only. ProbeST offers more flexibility to the probe-based chemistry by enabling the capture of both eukaryotic and prokaryotic transcript information.

ProbeST follows the [10X Genomics recommendations](https://www.10xgenomics.com/support/spatial-gene-expression-ffpe/documentation/workflows/ffpe-v-1/steps/experimental-design-and-planning/custom-probe-design-for-visium-spatial-gene-expression-and-chromium-single-cell-gene-expression-flex) for binding probes with the Visium platform.

<p align="center">
  <img width="700" alt="Chemistry_panel" src="https://github.com/user-attachments/assets/56b4686d-269a-4add-855c-b57d4f4e0861" />
</p>


ProbeST is built as a SnakeMake workflow, to be run from the terminal command line.
It is a modular workflow, as illustrated in the figure. The user can choose in the configuration file (config/config.yaml) different settings depending on the study.

## Case A: Probe design for a eukaryotic species
Choose pipeline version “1BLAST”. 
The optional rRNA filtering step is strongly recommended. 
For poorly annotated eukaryotic genomes, the optional cross-hybridization step is also strongly recommended.

## Case B: Host-pathogen studies 
Probe design for a prokaryotic species that you wish to detect in its host tissue (including host-pathogen interactions). 
Choose pipeline version “2BLASTS” and provide the CDS FASTA files of both the host genome and the prokaryotic genome. 
The rRNA filtering step is strongly recommended. 
Do not use the cross- hybridization step, as your probes are likely to be all filtered out by the databases used.

# Installation

First create a new directory with the name of your choice for Snakemake

```
mkdir directory_name
```

1. Install mamba

```
conda install -n base -c conda-forge mamba
```

2. Install Snakemake into an isolated software environment

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

3. Activate your Snakemake environment

```
conda activate snakemake
snakemake --help
snakemake --version
```


4. Install primer3

```
conda install -c bioconda primer3
```

5. Install blast

```
mamba install blast
```

6. Choose the configuration of choice and download **all** files in the 'workflow/' folder

Make sure the input CDS files are renamed so the Snakefile can recognize them:
-	pipeline version “1BLAST”: rename to CDS_gene_targets.fa the file with the genes of interest, and to CDS_all.fa the file with the full genome.
-	pipeline version “2BLASTS”: rename to CDS_gene_targets.fa the file with the genes of interest, to CDS_all_pathogen.fa and CDS_all_host.fa the full genomes of the pathogen and host genomes, respectively.
-	if choosing the rRNA filtering step: download the provided SILVA database. 
-	If choosing the cross-hybridisation step: 
  o	insert steps
  o	update the path in the cross_hyb.yaml configuration file to the path where your dbs are. 

   
8. Make a dry run of the workflow

```
snakemake -c 1 -s Snakefile_name --use-conda -np
```

The comment "missing output files" is normal. 


8. Run ProbeST workflow

```
snakemake -c 1 --use-conda -s Snakefile_name
```
or
```
snakemake --use-conda --conda-frontend mamba -c1
```

You can write
```
time snakemake -c 1 --use-conda -s Snakefile_name
```
if you want the total compute time.


9. If you want to delete outputs before running it again

```
snakemake -s Snakefile_name --delete-all-output
```

Check the log files in the logs/ folder, or the general Snakemake log file, for troubleshooting if an error occurs.


To build a DAG display (the pipeline workflow with the different steps/rules). You need to run ProbeST before in order to have both the inputs and outputs.

```
conda install graphviz
snakemake --dag selected_probes.txt | dot -Tsvg > dag.svg
```


## Input file formatting
In order to use this pipeline, the headers in the CDS_gene_targets.fa FASTA file with the genes of interest should be formatted accordingly: 

>key:value key:value key:value

Example:
>gene_ID:CBW18961.1 gene_name:sipC transcript_ID:CBW18961.1_2869

We provide an example of a python script to format the headers.


## For further information on Snakemake

Installation guide: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

Snakemake general tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html



# Probe design recommendations
ProbeST is aimed at designing custom probes for the Visium Spatial Gene Expression FFPE and Visium CytAssist Spatial Gene Expression assays.
Therefore, it follows the [current probe requirements](https://www.10xgenomics.com/support/spatial-gene-expression-ffpe/documentation/workflows/ffpe-v-1/steps/experimental-design-and-planning/custom-probe-design-for-visium-spatial-gene-expression-and-chromium-single-cell-gene-expression-flex) to optimise probe performance and specificity.
- A binding probe is 50 bp long, composed of two 25 bp sides: Left-Hand Side (LHS) and Right-Hand Side (RHS).
- Probes should be designed for coding regions of mRNA, rather than untranslated regions.
- GC content of each probe side should be between 44-72%.
- The 25th nucleotide (the 3' most nucleotide of the LHS) must be a T.
- Homopolymer repeats should be avoided.
- Off-target hybridisation should be checked for by BLASTing the probes against the reference transcriptome. This occurs if the probes bind to sequences other than the target mRNA sequence. Matches to off-target genes should have at least 5 mismatches in at least one of the LHS or RHS sides to prevent efficient hybridisation. Probes showing off-target hybridisation with less than or equal to 5 mismatches in either the LHS or the RHS should be removed.
- 3 probe pairs per target mRNA is recommended. If the gene is not long enough for 3 pairs, or if there are not 3 pairs specific enough, fewer probe pairs can be designed.
- There should not be any overlap between probe pair sequences to avoid competition between probes for the same binding site in the target mRNA.

# Cite it 
If you use this pipeline and the data, please cite the pre-print:  https://doi.org/10.1101/2025.08.29.673037

  
