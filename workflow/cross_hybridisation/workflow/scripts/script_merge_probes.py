# Open the file containing probe sequences
with open('probes.fasta', 'r') as file:
    lines = file.readlines()

# Check if the number of lines is a multiple of 4
if len(lines) % 4 != 0:
    exit(1)

# Open a new file to write the merged sequences
with open('merged_sequences.txt', 'w') as merged_file:
    # Process the lines in blocks of four (name, LHS, name, RHS)
    for i in range(0, len(lines), 4):
        lhs_name = lines[i].strip()  # Read the name for LHS
        lhs_sequence = lines[i + 1].strip()  # Read the LHS sequence
        rhs_name = lines[i + 2].strip()  # Read the name for RHS
        rhs_sequence = lines[i + 3].strip()  # Read the RHS sequence

        # Ensure that LHS and RHS sequences are from the same probe and have the same suffix
        lhs_parts = lhs_name.split('_LHS_')
        rhs_parts = rhs_name.split('_RHS_')

        if len(lhs_parts) == 2 and len(rhs_parts) == 2:
            lhs_probe_id, lhs_suffix = lhs_parts
            rhs_probe_id, rhs_suffix = rhs_parts

            if lhs_probe_id == rhs_probe_id and lhs_suffix == rhs_suffix:
                merged_sequence = lhs_sequence + rhs_sequence
                # Write the merged sequence to the file in FASTA format
                merged_file.write(f">{lhs_probe_id}_{lhs_suffix}\n{merged_sequence}\n")


