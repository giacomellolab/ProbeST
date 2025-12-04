from Bio import SeqIO
import pandas as pd
import sys
import generate_probe_pairs
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt


def load_fasta_probe_ids(fasta_file):
    """
    Extract probe identifiers from the FASTA file.

    Expected FASTA ID format:
        gene_LHS_<LHSID>_idx
        gene_RHS_<RHSID>_idx
    """
    fasta_entries = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        parts = record.id.split("_")
        gene_id = parts[0]
        side     = parts[3]  # LHS or RHS
        probe_id = parts[4]

        if gene_id not in fasta_entries:
            fasta_entries[gene_id] = {"LHS": set(), "RHS": set()}

        fasta_entries[gene_id][side].add(probe_id)

    return fasta_entries



def main(primer3_output_file, fasta_file, probes_csv, selected_probes, probe_quantifications, log_process):

    log = open(log_process, "w")


    # 1. Load original probe objects from primer3 output
    log.write("Reading primer3 output...\n")
    sequences = generate_probe_pairs.process(primer3_output_file)
    log.write(f"Loaded {len(sequences)} sequences from primer3\n")

    # 2. Load FASTA with filtered probes (that do not have significant hits to the db)
    log.write("Loading filtered FASTA probes...\n")
    fasta_probe_ids = load_fasta_probe_ids(fasta_file)
    log.write(f"Loaded filtered probes for {len(fasta_probe_ids)} genes\n")

    # 3. Keep only probes whose LHS and RHS are present in filtered fasta file
    log.write("Filtering probes based on FASTA entries...\n")

    for sequence in sequences:
        seq_gene = sequence.ID.split("_")[0]

        if seq_gene not in fasta_probe_ids:
            # No probes retained for this gene
            sequence.PROBES = []
            continue

        keep_LHS = fasta_probe_ids[seq_gene]["LHS"]
        keep_RHS = fasta_probe_ids[seq_gene]["RHS"]

        filtered = []
        for probe in sequence.PROBES:
            lhs_id = probe.LHS_ID
            rhs_id = probe.RHS_ID

            if (lhs_id in keep_LHS) and (rhs_id in keep_RHS):
                filtered.append(probe)

        sequence.PROBES = filtered
        log.write(f"{sequence.ID}: {len(sequence.PROBES)} probes retained after FASTA filter\n")

  
    # 4. Remove RHS homopolymers > 5
    log.write("Applying homopolymer filter...\n")

    MAX_HOMOPOLYMER = 5
    for sequence in sequences:
        filtered = []
        for probe in sequence.PROBES:
            rhs = probe.RHS[:25]
            remove = False
            for base in "ATGC":
                if base * MAX_HOMOPOLYMER in rhs:
                    remove = True
                    break
            if not remove:
                filtered.append(probe)

        sequence.PROBES = filtered
        log.write(f"{sequence.ID}: {len(sequence.PROBES)} after homopolymer filter\n")


    # 5. Remove overlapping probes
    log.write("Applying overlap filter...\n")

    for sequence in sequences:
        i = 0
        while i < len(sequence.PROBES) - 1:
            if (int(sequence.PROBES[i].START) + 4) > int(sequence.PROBES[i+1].END):
                del sequence.PROBES[i+1]
            else:
                i += 1
        log.write(f"{sequence.ID}: {len(sequence.PROBES)} after overlap filter\n")

    # 6. Keep first 3 probe pairs
    log.write("Limiting to first 3 probes per sequence...\n")

    for sequence in sequences:
        sequence.PROBES = sequence.PROBES[:3]
        log.write(f"{sequence.ID}: retained {len(sequence.PROBES)} probes\n")


    # 7. Generate probe CSV (SpaceRanger input)
    log.write("Generating probe CSV...\n")

    probe_ref = pd.DataFrame(columns=["gene_id", "probe_seq", "probe_id", "included"])

    for sequence in sequences:
        seq_id, gene_name = sequence.ID.split("_")[0], sequence.ID.split("_")[1]

        for probe in sequence.PROBES:
            lhs = probe.LHS[-25:]
            rhs = probe.RHS[:25]
            full_probe = lhs + rhs
            probe_id = f"{seq_id}|{gene_name}|{probe.HASH_ID}"

            df_temp = pd.DataFrame(
                {
                    "gene_id": seq_id,
                    "probe_seq": full_probe,
                    "probe_id": probe_id,
                    "included": "TRUE",
                },
                index=[0],
            )

            probe_ref = pd.concat([probe_ref, df_temp], ignore_index=True)

    probe_ref = probe_ref.drop_duplicates().sort_values("probe_id")
    probe_ref.to_csv(probes_csv, index=False)

    # 8. Save selected probes txt file
    file = open(selected_probes, "w")
    sys.stdout = file
    generate_probe_pairs.printer(sequences)
    file.close()

    # 9. Probe quantification report
    log.write("Generating probe quantification report...\n")
    def probe_counter(sequences):
        n = 0
        for sequence in sequences: 
            for probe in sequence.PROBES:
                n = n + 1 
        return n

    # Count sequences by number of probes
    count_occurrences = {} 
    for sequence in sequences:
        sequence.PROBES_count = len(sequence.PROBES)
        if sequence.PROBES_count in count_occurrences:
            count_occurrences[sequence.PROBES_count] += 1
        else:
            count_occurrences[sequence.PROBES_count] = 1

    count_occ = {}
    for s in sequences:
        c = len(s.PROBES)
        count_occ[c] = count_occ.get(c, 0) + 1

    # Make probe quantifications file
    file = open(probe_quantifications, "w")
    sys.stdout = file

    print("Total probes:", probe_counter(sequences), "\n")
    print("the total amount of sequences designed probes for: ", len(sequences)- count_occurrences[0])
    for i in range(0,5): #adjust range of to the maximum amount of probes per gene which you are interested in
        if i in count_occurrences:
            print("the amount of sequences containing ", i, " probes: ", count_occurrences[i])    
        else:
            print('the amount of sequences containing ', i, ' probes:  0')

    # melting temperature analysis
    print("\nTm consistency analysis:\n")
    tm_per_gene = {}

    for sequence in sequences:
        seq_id, gene_name = sequence.ID.split("_")[0], sequence.ID.split("_")[1]

        if gene_name not in tm_per_gene:
            tm_per_gene[gene_name] = []

        for probe in sequence.PROBES:
            lhs = probe.LHS[-25:]
            rhs = probe.RHS[:25]
            full_probe = lhs + rhs
            tm = mt.Tm_NN(full_probe)
            tm_per_gene[gene_name].append(tm)

    tm_ranges = {gene: (max(tms) - min(tms)) if len(tms) > 1 else 0 for gene, tms in tm_per_gene.items()}

    all_ranges = list(tm_ranges.values())
    print(f"Average ΔTm per gene: {np.mean(all_ranges):.2f} °C")
    print(f"Max ΔTm: {np.max(all_ranges):.2f} °C\n")

    flagged = []
    for gene, delta in tm_ranges.items():
        tms = tm_per_gene[gene]
        print(f"Gene {gene}")
        print(f"  Tms: {', '.join(f'{x:.2f}' for x in tms)}")
        print(f"  ΔTm = {delta:.2f} °C")
        if delta > 5:
            print("  WARNING: ΔTm > 5 °C\n")
            flagged.append((gene, delta))
        print()

    if flagged:
        print("\nGenes with ΔTm > 5 °C:")
        for g, d in flagged:
            print(f"{g}: {d:.2f} °C")

    file.close()

    log.close()

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python3 select_probe_pairs.py <primer3_output_file> <fasta_file> <probes_csv> <selected_probes> <probe_quantifications> <log_process> ")
        sys.exit(1)

    primer3_output_file = sys.argv[1]
    fasta_file = sys.argv[2]
    probes_csv = sys.argv[3]
    selected_probes = sys.argv[4]
    probe_quantifications = sys.argv[5]
    log_process = sys.argv[6]

    main(primer3_output_file, fasta_file, probes_csv, selected_probes, probe_quantifications, log_process)