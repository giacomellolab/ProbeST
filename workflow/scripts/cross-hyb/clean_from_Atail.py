import sys

# Function to trim poly-A tail from a sequence
def trim_poly_a_tail(sequence):
    return sequence.rstrip('A')

def main(input_file, output_file):
    # Open the merged sequences file
    with open(input_file, 'r') as merged_file:
        lines = merged_file.readlines()

    # Open a new file to write the cleaned sequences
    with open(output_file, 'w') as cleaned_file:
        # Process the lines in pairs (name, sequence)
        for i in range(0, len(lines), 2):
            name = lines[i].strip()  # Read the name line
            sequence = lines[i + 1].strip()  # Read the sequence line

            # Trim the poly-A tail from the sequence
            cleaned_sequence = trim_poly_a_tail(sequence)

            # Write the cleaned sequence to the new file in FASTA format
            cleaned_file.write(f"{name}\n{cleaned_sequence}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise ValueError("Usage: python3 clean_from_Atail.py input.txt output.fa")

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    main(input_file, output_file)
