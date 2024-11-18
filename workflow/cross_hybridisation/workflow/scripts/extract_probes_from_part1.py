def extract_probes_from_file(file_path):
    probes = []
    with open(file_path, 'r') as file:
        probe_info = {}
        for line in file:
            line = line.strip()
            if line.startswith("probe pair LHS ID:"):
                if probe_info:
                    probes.append(probe_info)
                probe_info = {}
                probe_info['LHS ID'] = line.split(": ")[1]
            elif line.startswith("probe pair RHS ID:"):
                probe_info['RHS ID'] = line.split(": ")[1]
            elif line.startswith("probe pair START:"):
                probe_info['START'] = line.split(": ")[1]
            elif line.startswith("probe pair END:"):
                probe_info['END'] = line.split(": ")[1]
            elif line.startswith("probe LHS:"):
                probe_info['LHS sequence'] = line.split(": ")[1]
            elif line.startswith("probe RHS:"):
                probe_info['RHS sequence'] = line.split(": ")[1]
        if probe_info:
            probes.append(probe_info)
    return probes

def write_fasta(probes, output_file):
    with open(output_file, 'w') as file:
        for probe in probes:
            if "LHS sequence" in probe:
                file.write(f'>{probe["LHS ID"]}\n{probe["LHS sequence"]}\n')
            if "RHS sequence" in probe:
                file.write(f'>{probe["RHS ID"]}\n{probe["RHS sequence"]}\n')

if __name__ == "__main__":
    probes = extract_probes_from_file("selected_probes.txt")
    write_fasta(probes, "probes.fasta")
    #print("FASTA file generated successfully!")


