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
    "CH": {"C": 1, "H": 0},
    "carbonyl": {"C": 1, "O": 1},
    "amino": {"N": 1, "H": 1},
}

def format_structure(component):
    #Formats the CHNOP count of a component into a string.
    return "".join(f"{element}{count}" for element, count in structures[component].items())

def calculate_full_peptide(sequence):
    total_formula = Counter()
    format_formula = []

    total_formula.update(structures["N_term"])
    format_formula.append(format_structure("N_term"))
    total_formula.update(structures["CH"])
    format_formula.append(format_structure("CH"))

    for residue in sequence:
        if residue not in structures:
            raise ValueError(f"Invalid residue: {residue}")
        
        total_formula.update(structures[residue])
        format_formula.append(format_structure(residue))
        total_formula.update(structures["carbonyl"])
        format_formula.append(format_structure("carbonyl"))
        total_formula.update(structures["amino"])
        format_formula.append(format_structure("amino"))

    total_formula.update(structures["C_term"])
    format_formula.append(format_structure("C_term"))
    return format_formula, total_formula

def simulate_fragmentation(sequence):
    fragments = {"a": [], "b": [], "c": [], "x": [], "y": [], "z": []}
    fragment_structures = {"a": [], "b": [], "c": [], "x": [], "y": [], "z": []}

    for i, residue in enumerate(sequence):
        
        if residue not in structures:
            raise ValueError(f"Invalid nucleotide: {residue}")
        
        # N-terminal side 
        prefix_formula = Counter()
        prefix_structure = []

        prefix_formula.update(structures["N_term"])
        prefix_structure.append(format_structure("N_term"))

        for j in range(i + 1):
            prefix_structure.append(format_structure("CH"))
            prefix_formula.update(structures["CH"])
            prefix_structure.append(format_structure(sequence[j]))
            prefix_formula.update(structures[sequence[j]])
            prefix_structure.append(format_structure("carbonyl"))
            prefix_formula.update(structures["carbonyl"])
            prefix_structure.append(format_structure("amino"))
            prefix_formula.update(structures["amino"])
        
        fragments["a"].append(prefix_formula - Counter(structures["carbonyl"]) - Counter(structures["amino"]))
        fragment_structures["a"].append("--".join(prefix_structure + [format_structure("carbonyl")] + [format_structure("amino")]))

        fragments["b"].append(prefix_formula - Counter(structures["amino"]))
        fragment_structures["b"].append("--".join(prefix_structure + [format_structure("amino")]))

        fragments["c"].append(prefix_formula)
        fragment_structures["c"].append("--".join(prefix_structure))
        
        if i == len(sequence)-1:
            break
        
        # C-terminal side
        suffix_formula = Counter()
        suffix_structure = []

        for j in range(i+1, len(sequence)):
            suffix_structure.append(format_structure("carbonyl"))
            suffix_formula.update(structures["carbonyl"])
            suffix_structure.append(format_structure("amino"))
            suffix_formula.update(structures["amino"])
            suffix_structure.append(format_structure("CH"))
            suffix_formula.update(structures["CH"])
            suffix_structure.append(format_structure(sequence[j]))
            suffix_formula.update(structures[sequence[j]])

            if j == len(sequence)-1:
                suffix_structure.append(format_structure("C_term"))
                suffix_formula.update(structures["C_term"])

        fragments["z"].append(suffix_formula - Counter(structures["carbonyl"]) - Counter(structures["amino"]))
        fragment_structures["z"].append("--".join(suffix_structure + [format_structure("carbonyl")] + [format_structure("amino")]))

        fragments["y"].append(suffix_formula - Counter(structures["carbonyl"]))
        fragment_structures["y"].append("--".join(suffix_structure + [format_structure("carbonyl")]))

        fragments["x"].append(suffix_formula)
        fragment_structures["x"].append("--".join(suffix_structure))

    return fragments, fragment_structures

def format_chemical_formula(total_formula):
    return ''.join(f"{element}{count}" for element, count in sorted(total_formula.items()))

def main():
    parser = argparse.ArgumentParser(description="Simulate peptide fragmentation into b, y, c, x, and z ions.")
    parser.add_argument("sequence", type=str, help="Peptide sequence (e.g., ACDEFGHIK)")
    args = parser.parse_args()

    sequence = args.sequence.upper()

    try:
        # Calculate the full peptide formula
        formatted_structure, total_formula = calculate_full_peptide(sequence)
        chemical_formula = format_chemical_formula(total_formula)

        print(f"Sequence to fragment: {sequence}")
        print(f"\ntotal_formula: {total_formula}")
        print(f"\nChemical Formula: {chemical_formula}")
        print(f"\n{sequence} Linkage Structure:")
        print("--".join(formatted_structure))

        # Simulate fragmentation
        fragments, fragment_structures = simulate_fragmentation(sequence)

        for ion_type, ion_list in fragments.items():
            print(f"\n{ion_type}-series fragments:")
            if ion_type in ["w", "x", "y", "z"]:
                # Reverse numbering for w, x, y, z series
                total_fragments = len(ion_list)
                for i, fragment in enumerate(ion_list):
                    fragment_number = total_fragments - i  # Reverse numbering
                    formatted_fragment = format_chemical_formula(fragment)
                    print(f"Fragment {fragment_number} Linked Structure:")
                    print(fragment_structures[ion_type][i])
                    print(f"Fragment {fragment_number} Formula: {formatted_fragment}")

            else:
                # Keep original numbering for a, b, c, d series
                for i, fragment in enumerate(ion_list):
                    fragment_number = i + 1  # Regular numbering
                    formatted_fragment = format_chemical_formula(fragment)
                    print(f"Fragment {fragment_number} Linked Structure:")
                    print(fragment_structures[ion_type][i])
                    print(f"Fragment {fragment_number} Formula: \n{formatted_fragment}")

    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()

