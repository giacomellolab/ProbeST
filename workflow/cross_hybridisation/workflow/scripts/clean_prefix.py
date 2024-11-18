# Function to trim the specified prefix from a sequence
def trim_prefix(sequence, prefix):
    if sequence.startswith(prefix):
        return sequence[len(prefix):]
    return sequence

# The prefix sequence to be trimmed
prefix_sequence = 'CCTTGGCACCCGAGAATTCCA'

# Open the cleaned sequences file (input)
with open('cleaned_sequences.fasta', 'r') as cleaned_input_file:
    lines = cleaned_input_file.readlines()

# Open a new file to write the sequences cleaned from the prefix (output)
with open('cleaned_sequences_from_prefix.fasta', 'w') as cleaned_output_file:
    # Process the lines in pairs (name, sequence)
    for i in range(0, len(lines), 2):
        name = lines[i].strip().replace(">>", ">")  # Read the name line and remove one >
        sequence = lines[i + 1].strip()  # Read the sequence line

        # Trim the prefix from the sequence
        cleaned_sequence = trim_prefix(sequence, prefix_sequence)

        # Write the cleaned sequence to the new file in FASTA format
        cleaned_output_file.write(f"{name}\n{cleaned_sequence}\n")

