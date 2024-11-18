# Function to trim poly-A tail from a sequence
def trim_poly_a_tail(sequence):
    return sequence.rstrip('A')

# Open the merged sequences file
with open('merged_sequences.txt', 'r') as merged_file:
    lines = merged_file.readlines()

# Open a new file to write the cleaned sequences
with open('cleaned_sequences.fasta', 'w') as cleaned_file:
    # Process the lines in pairs (name, sequence)
    for i in range(0, len(lines), 2):
        name = lines[i].strip()  # Read the name line
        sequence = lines[i + 1].strip()  # Read the sequence line

        # Trim the poly-A tail from the sequence
        cleaned_sequence = trim_poly_a_tail(sequence)

        # Write the cleaned sequence to the new file in FASTA format
        cleaned_file.write(f"{name}\n{cleaned_sequence}\n")

# Notify user of completion
#print("Cleaning completed. Cleaned sequences are saved in 'cleaned_sequences.fasta'.")

