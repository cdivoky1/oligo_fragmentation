import argparse
from collections import Counter

# Define the structures of amino acids and relevant components
structures = {
    "A": {"C": 0, "H": 1, "N": 0, "O": 0},  # Alanine ALA A
    "R": {"C": 6, "H": 12, "N": 4, "O": 1}, # Arginine
    "N": {"C": 4, "H": 6, "N": 2, "O": 2},  # Asparagine
    "D": {"C": 4, "H": 5, "N": 1, "O": 3},  # Aspartic acid
    "C": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1}, # Cysteine
    "E": {"C": 5, "H": 7, "N": 1, "O": 3},  # Glutamic acid
    "Q": {"C": 5, "H": 8, "N": 2, "O": 2},  # Glutamine
    "G": {"C": 2, "H": 3, "N": 1, "O": 1},  # Glycine
    "H": {"C": 6, "H": 7, "N": 3, "O": 1},  # Histidine
    "I": {"C": 6, "H": 11, "N": 1, "O": 1}, # Isoleucine
    "L": {"C": 6, "H": 11, "N": 1, "O": 1}, # Leucine
    "K": {"C": 6, "H": 12, "N": 2, "O": 1}, # Lysine
    "M": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1}, # Methionine
    "F": {"C": 9, "H": 9, "N": 1, "O": 1},  # Phenylalanine
    "P": {"C": 5, "H": 7, "N": 1, "O": 1},  # Proline
    "S": {"C": 3, "H": 5, "N": 1, "O": 2},  # Serine
    "T": {"C": 4, "H": 7, "N": 1, "O": 2},  # Threonine
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},# Tryptophan
    "Y": {"C": 9, "H": 9, "N": 1, "O": 2},  # Tyrosine
    "V": {"C": 5, "H": 9, "N": 1, "O": 1},  # Valine
    "N_term": {"C": 0, "H": 2, "N": 1, "O": 0},  # N_terminal residue
    "C_term": {"C": 1, "H": 1, "N": 0, "O": 2},  # C_terminal residue
    "H": {"H": 1},
    "C": {"C": 1},
}

def format_chemical_formula(counter):
    return "".join(f"{element}{count}" for element, count in sorted(counter.items()))

def calculate_full_peptide(sequence):
    total_formula = Counter(structures["N_term"])

    for residue in sequence:
        if residue not in structures:
            raise ValueError(f"Invalid residue: {residue}")
        total_formula.update(structures[residue])

    total_formula.update(structures["C_term"])
    return total_formula

def simulate_fragmentation(sequence):
    fragments = {"b": [], "y": [], "c": [], "x": [], "z": []}

    for i in range(len(sequence)):
        # Generate b and c fragments (N-terminal)
        b_formula = Counter(structures["N_term"])
        c_formula = Counter(structures["N_term"])

        for residue in sequence[:i + 1]:
            b_formula.update(structures[residue])
            c_formula.update(structures[residue])

        c_formula.update({"H": 1, "O": 1})  # Add H and O for c-fragments
        fragments["b"].append(format_chemical_formula(b_formula))
        fragments["c"].append(format_chemical_formula(c_formula))

        # Generate y, x, z fragments (C-terminal)
        y_formula = Counter(structures["C_term"])
        x_formula = Counter(structures["C_term"])
        z_formula = Counter(structures["C_term"])

        for residue in sequence[i:]:
            y_formula.update(structures[residue])
            x_formula.update(structures[residue])
            z_formula.update(structures[residue])

        x_formula.update({"O": 1})  # Add O for x-fragments
        z_formula.update({"H": -1})  # Remove H for z-fragments

        fragments["y"].append(format_chemical_formula(y_formula))
        fragments["x"].append(format_chemical_formula(x_formula))
        fragments["z"].append(format_chemical_formula(z_formula))

    return fragments

def main():
    parser = argparse.ArgumentParser(description="Simulate peptide fragmentation into b, y, c, x, and z ions.")
    parser.add_argument("sequence", type=str, help="Peptide sequence (e.g., ACDEFGHIK)")
    args = parser.parse_args()

    sequence = args.sequence.upper()

    try:
        # Calculate the full peptide formula
        full_formula = calculate_full_peptide(sequence)
        print("Full Peptide Formula:", format_chemical_formula(full_formula))

        # Simulate fragmentation
        fragments = simulate_fragmentation(sequence)

        for ion_type, ion_list in fragments.items():
            print(f"\n{ion_type}-series fragments:")
            for i, fragment in enumerate(ion_list, 1):
                print(f"Fragment {i}: {fragment}")

    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()

