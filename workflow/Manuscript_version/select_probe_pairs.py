"""
This is a script for selecting final probe pairs based on: 
- off target hybridisation (BLAST OUTPUT)
- max homopolymer repeats RHS = 5 
- max probe pairs per gene = 3

Input: 
- python code for generating initial probe pairs 
- trimmed BLAST output generated using target_specificity_trim.py 

Output: 
- final probe pairs overview txt file 
- probe quantifications txt file 
- probe reference csv file compatible with Space Ranger
"""

import pandas as pd 
import sys
import generate_probe_pairs_Snakemake_pipeline



def main(primer3_output_file, blast_output_file, probes_csv, selected_probes, probe_quantifications, log_process):
    # Open a log file to record debug information
    log_file = open(log_process, 'w')

    # Import initial probe pairs 
    log_file.write("Reading primer3 output file...\n")
    sequences = generate_probe_pairs_Snakemake_pipeline.process(primer3_output_file)
    log_file.write(f"Number of sequences processed: {len(sequences)}\n")

    # Import BLAST output for off-target hybridisation
    log_file.write("Reading BLAST output file...\n")
    hits = pd.read_csv(blast_output_file, sep = ' ')
    log_file.write(f"Number of BLAST hits: {len(hits)}\n")

    # Delete probes which have >1 BLAST hits on both the RHS and LHS (non-specific probes) Based on the 10X guidelines, at least one of the two probe sides needs to have >5 mismatches to prevent off-target hybridisation  
    hits['probe_side'] = hits['probe_id'].str.split("_").str[4]  #create column indicating RHS or LHS (index 4 for S.Tm, but 3 is default as: geneid_genename_transcriptid_probeside)
    #hits['probe_side'] = hits['probe_id'].str.split("_").str[3] #create column indicating RHS or LHS (index 4 for S.Tm, but 3 is default as: geneid_genename_transcriptid_probeside)
    grouped = hits.groupby(['probe_pair_id', 'probe_side']).size().reset_index(name='count')
    filtered_groups = grouped.pivot(index='probe_pair_id', columns='probe_side', values='count')
    filtered_groups = filtered_groups[(filtered_groups['LHS'] > 2) & (filtered_groups['RHS'] > 2)]
    nonspec_probes = filtered_groups.index.tolist() #create lists of probes which have off-targets for both RHS and LHS probes 
    log_file.write(f"Number of non-specific probes: {len(nonspec_probes)}\n")

    # # Alternative: Delete probes which have >1 BLAST hit (5 or less mismatches) on one or both probe sides
    # hits2 = hits[hits['probe_id'].duplicated(keep = False) == True] #make dataframe only containing probes with >1 hit
    # nonspec_probes = hits2['probe_id'].values.tolist() #create list of probe_ids which have >1 hits (non-specific) 

    for sequence in sequences: 
        for i, probe in enumerate(sequence.PROBES): 
                for nonspec_probe in nonspec_probes:
                    if nonspec_probe == probe.LHS_ID: 
                        del sequence.PROBES[i]
                        break 
        log_file.write(f"Number of probes after specificity filter for sequence {sequence.ID}: {len(sequence.PROBES)}\n")


    # Delete probes (RHS) which contain a homopolymer repeat of >5 identical nucleotides 
    K = 5 #initiate max length of homopolymer repeat 
    for sequence in sequences: 
        for i, probe in enumerate(sequence.PROBES):
            hyb_RHS = probe.RHS[:25]
            for idx, ele in enumerate(hyb_RHS):
                substr = ele * K #create equi string
                if hyb_RHS[idx : idx + K] > substr:
                    del sequence.PROBES[i]
                    break


    # Delete probes which overlap 
    for sequence in sequences:
        iterator = 0
        while iterator < len(sequence.PROBES)-1:
            if (int(sequence.PROBES[iterator].END) + 4)> int(sequence.PROBES[iterator+1].START): #have a minimum space of 4 nucleotides between two probe pairs (hyb regions). Somewhat arbitrarily chosen. 
                del sequence.PROBES[iterator+1]
            else:
                iterator += 1
        log_file.write(f"Number of probes after overlap filter for sequence {sequence.ID}: {len(sequence.PROBES)}\n")

    # Keep first 3 probe pairs 
    for index, sequence in enumerate(sequences): 
        while len(sequence.PROBES) > 3:
            del sequence.PROBES[-1]
        log_file.write(f"Number of probes after limiting to 3 for sequence {sequence.ID}: {len(sequence.PROBES)}\n")

    
    # Generate the 3 output files

    # Create SpaceRanger compatible CSV reference file.
    # "Region" Column not included because probes have not been designed specifically for bridging splice junctions 
    probe_ref = pd.DataFrame(columns = ["gene_id", "probe_seq", "probe_id", "included"])
    x = 1 
    for sequence in sequences: 
        #gene_name = "gene_" + str(x) #generate incrementing gene names, in case you don't have gene name annotations for your genome
        for probe in sequence.PROBES:
            hyb_LHS = probe.LHS[-25:] #only take the hybrdising part of the probe, i.e. the last 25 nucleotides    
            hyb_RHS = probe.RHS[:25] #only take the hybridising part of the probe, i.e. the first 25 nucleotides 
            full_probe = hyb_LHS + hyb_RHS 
            sequence_id = sequence.ID.split('_')[0]
            gene_name = sequence.ID.split('_')[1]
            probe_id = str(sequence_id + "|" + gene_name + "|" + probe.HASH_ID)
            df_temp = pd.DataFrame({"gene_id": sequence_id, "probe_seq": full_probe, "probe_id": probe_id, "included": "TRUE"}, index = [0])
            probe_ref = probe_ref.drop_duplicates()
            probe_ref = pd.concat([df_temp, probe_ref.loc[:]]).reset_index(drop=True)
        x = x + 1
    probe_ref = probe_ref.sort_values(by = ['probe_id'])
    probe_ref.to_csv(probes_csv, index = False) #export dataframe to csv

    # Save overview of generated probes in txt file 
    file = open(selected_probes, "w")
    sys.stdout = file
    generate_probe_pairs_Snakemake_pipeline.printer(sequences)
    file.close()


    # Count final amount of probes, for the last file
    def probe_counter(sequences):
        n = 0
        for sequence in sequences: 
            for probe in sequence.PROBES:
                n = n + 1 
        return n

    # Create dictionary counting the occurrences of n probes per gene, for the last file   
    count_occurrences = {} 
    for sequence in sequences:
        sequence.PROBES_count = len(sequence.PROBES)
        if sequence.PROBES_count in count_occurrences:
            count_occurrences[sequence.PROBES_count] += 1
        else:
            count_occurrences[sequence.PROBES_count] = 1

    # Save probe counts and other quantifications in txt file 
    file = open(probe_quantifications, "w")
    sys.stdout = file 
    print("the final amount of probes is: ", probe_counter(sequences), "\n")   
    print("the total amount of sequences designed probes for: ", len(sequences)- count_occurrences[0])
    for i in range(0,5): #adjust range of to the maximum amount of probes per gene which you are interested in
        if i in count_occurrences:
            print("the amount of sequences containing ", i, " probes: ", count_occurrences[i])    
        else:
            print('the amount of sequences containing ', i, ' probes:  0')
    file.close()


    # Close the log file
    log_file.close()

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python3 select_probe_pairs_Snakemake_pipeline.py <primer3_output_file> <blast_output_file> <probes_csv> <selected_probes> <probe_quantifications> <log_process> ")
        sys.exit(1)

    primer3_output_file = sys.argv[1]
    blast_output_file = sys.argv[2]
    probes_csv = sys.argv[3]
    selected_probes = sys.argv[4]
    probe_quantifications = sys.argv[5]
    log_process = sys.argv[6]

    main(primer3_output_file, blast_output_file, probes_csv, selected_probes, probe_quantifications, log_process)


