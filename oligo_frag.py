import argparse
from collections import Counter

# Define the structures
structures = {
    "sugar": {"C": 5, "H": 9, "O": 2}, #3' OH intentially left off
    "outer_link": {"O": 1},
    "inner_link": {"H": 1, "O": 2, "P": 1},
    "A": {"C": 5, "H": 4, "N": 5, "O": 0},
    "U": {"C": 4, "H": 4, "N": 2, "O": 2},
    "C": {"C": 4, "H": 4, "N": 3, "O": 1},
    "G": {"C": 5, "H": 4, "N": 3, "O": 1},
    "H": {"H": 1},  # For 3'-H
    "OH": {"O":1, "H": 1},  # For 3'-OH
}

def format_structure(component):
    """Formats the CHNOP count of a component into a string."""
    return "".join(f"{element}{count}" for element, count in structures[component].items())

def calculate_chemical_formula(sequence):
    total_formula = Counter()
    formatted_structure = []

    for nucleotide in sequence:
        if nucleotide not in structures:
            raise ValueError(f"Invalid nucleotide: {nucleotide}")

        # Add the components for this nucleotide
        formatted_structure.append(format_structure("outer_link"))
        total_formula.update(structures["outer_link"])

        formatted_structure.append(format_structure("inner_link"))
        total_formula.update(structures["inner_link"])

        formatted_structure.append(format_structure("outer_link"))
        total_formula.update(structures["outer_link"])

        formatted_structure.append(format_structure("sugar"))
        total_formula.update(structures["sugar"])

        formatted_structure.append(format_structure(nucleotide))
        total_formula.update(structures[nucleotide])

    # Handle final outer_link
    formatted_structure.append(format_structure("outer_link"))
    total_formula.update(structures["outer_link"])

    # Account for 3'OH on terminal end
    formatted_structure.append(format_structure("H"))
    total_formula.update(structures["H"])

    return formatted_structure, total_formula

def simulate_fragmentation(sequence):
    """Simulates RNA fragmentation into a, b, c, d, w, x, y, z ions."""

    fragments = {"a": [], "a-bh": [], "b": [], "c": [], "d": [], "w": [], "x": [], "y": [], "z": []}
    fragment_structures = {"a": [], "a-bh": [], "b": [], "c": [], "d": [], "w": [], "x": [], "y": [], "z": []}

    for i, nucleotide in enumerate(sequence):

        if nucleotide not in structures:
            raise ValueError(f"Invalid nucleotide: {nucleotide}")

        # 5' side fragments (a, b, c, d)
        prefix_formula = Counter()
        prefix_structure = []
        for j in range(i + 1):
            prefix_structure.append(format_structure("outer_link"))
            prefix_formula.update(structures["outer_link"])
            prefix_structure.append(format_structure("inner_link"))
            prefix_formula.update(structures["inner_link"])
            prefix_structure.append(format_structure("outer_link"))
            prefix_formula.update(structures["outer_link"])
            prefix_structure.append(format_structure("sugar"))
            prefix_formula.update(structures["sugar"])
            prefix_structure.append(format_structure(sequence[j]))
            prefix_formula.update(structures[sequence[j]])

        # Add 3'H group for 'a' fragments
        fragments["a"].append(prefix_formula + Counter(structures["H"]))
        fragment_structures["a"].append("--".join(prefix_structure + [format_structure("H")])) 

        # Remove base (last structure) and add a 3' H onto sugar
        fragments["a-bh"].append(prefix_formula - Counter(structures[sequence[i]]) + Counter(structures["H"]))
        fragment_structures["a-bh"].append("--".join(prefix_structure[:-1] + [format_structure("H")])) 

        # Add 3' OH
        fragments["b"].append(prefix_formula + Counter(structures["OH"]))
        fragment_structures["b"].append("--".join(prefix_structure + [format_structure("OH")]))

        # Add 3' O and PO2H (unsure whether the phoshate is protonated once more in the machine during fragmentation)
        fragments["c"].append(prefix_formula + Counter(structures["outer_link"]) + Counter(structures["inner_link"]))
        fragment_structures["c"].append("--".join(prefix_structure + [format_structure("outer_link")] + [format_structure("inner_link")]))

        # Add 3' Phosphate group (PO4) with an added hydogen on the terminal O of the PO3 group
        fragments["d"].append(prefix_formula + Counter(structures["outer_link"]) + Counter(structures["inner_link"]) + Counter(structures["outer_link"]) + Counter(structures["H"]))
        fragment_structures["d"].append("--".join(prefix_structure + [format_structure("outer_link")] + [format_structure("inner_link")] + [format_structure("outer_link")] + [format_structure("H")]))

        # 3' side fragments (w, x, y, z)
        suffix_formula = Counter()
        suffix_structure = []
        for j in range(i+1, len(sequence)): # j starts at i+1 so it doesn't capture preceding base
            suffix_structure.append(format_structure(sequence[j]))
            suffix_formula.update(structures[sequence[j]])
            suffix_structure.append(format_structure("sugar"))
            suffix_formula.update(structures["sugar"])
            if j < len(sequence)-1:
                suffix_structure.append(format_structure("outer_link"))
                suffix_formula.update(structures["outer_link"])
                suffix_structure.append(format_structure("inner_link"))
                suffix_formula.update(structures["inner_link"])
                suffix_structure.append(format_structure("outer_link"))
                suffix_formula.update(structures["outer_link"])
            if j == len(sequence)-1:
                suffix_structure.append(format_structure("OH"))
                suffix_formula.update(structures["OH"])
            
        fragments["w"].append(suffix_formula - Counter(structures["sugar"]))
        fragment_structures["w"].append("--".join(suffix_structure))  

        fragments["x"].append(suffix_formula)
        fragment_structures["x"].append("--".join(suffix_structure))

        fragments["y"].append(suffix_formula + Counter(structures["inner_link"]))
        fragment_structures["y"].append("--".join(suffix_structure + [format_structure("inner_link")]))

        fragments["z"].append(suffix_formula + Counter(structures["outer_link"]))
        fragment_structures["z"].append("--".join(suffix_structure + [format_structure("outer_link")]))

    return fragments, fragment_structures

def format_chemical_formula(total_formula):
    return ''.join(f"{element}{count}" for element, count in sorted(total_formula.items()))

def main():
    parser = argparse.ArgumentParser(description="Build RNA oligo and simulate fragmentation.")
    parser.add_argument("sequence", type=str, help="RNA sequence (e.g., ACGU)")
    parser.add_argument("--output", type=str, default="chemical_formula.txt", help="Output file for chemical formula")
    args = parser.parse_args()

    sequence = args.sequence.upper()
    try:
        # Build full RNA molecule
        formatted_structure, total_formula = calculate_chemical_formula(sequence)
        print("Linkage Structure:")
        print("--".join(formatted_structure))

        chemical_formula = format_chemical_formula(total_formula)
        print(f"total_formula: \n{total_formula}")
        print("\nChemical Formula:")
        print(chemical_formula)

        # Simulate fragmentation
        fragments, fragment_structures = simulate_fragmentation(sequence)
        print("\nFragmentation Results:")
        for ion_type, ion_list in fragments.items():
            print(f"\n{ion_type}-series fragments:")
            for i, fragment in enumerate(ion_list):
                formatted_fragment = format_chemical_formula(fragment)
                print(f"Fragment {i + 1} Linked Structure:")
                print(fragment_structures[ion_type][i])
                print(f"Fragment {i + 1} Formula:")
                print(formatted_fragment)

        with open(args.output, "w") as file:
            file.write(chemical_formula + "\n")
        print(f"\nChemical formula saved to {args.output}")
    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()
