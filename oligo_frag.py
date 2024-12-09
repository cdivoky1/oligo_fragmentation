import argparse
from collections import Counter

# Define the structures
structures = {
    "sugar": {"C": 5, "H": 9, "O": 2},
    "outer_link": {"O": 1},
    "inner_link": {"H": 1, "O": 2, "P": 1},
    "A": {"C": 5, "H": 4, "N": 5, "O": 0},
    "U": {"C": 4, "H": 4, "N": 2, "O": 2},
    "C": {"C": 4, "H": 4, "N": 3, "O": 1},
    "G": {"C": 5, "H": 4, "N": 3, "O": 1},
    "H": {"H": 1},  # For 3'-OH
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
    fragments = {"a": [], "b": [], "c": [], "d": [], "w": [], "x": [], "y": [], "z": []}
    fragment_structures = {"a": [], "b": [], "c": [], "d": [], "w": [], "x": [], "y": [], "z": []}

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
            prefix_structure.append(format_structure(sequence[i]))
            prefix_formula.update(structures[sequence[i]])
            #if j < i:  # Add inner/outer linkers for all but the last nucleotide
            #    prefix_structure.append(format_structure("outer_link"))
            #    prefix_formula.update(structures["outer_link"])
            #    prefix_structure.append(format_structure("inner_link"))
            #    prefix_formula.update(structures["inner_link"])
        fragments["a"].append(prefix_formula - Counter(structures["sugar"]))
        fragment_structures["a"].append("--".join(prefix_structure[:-1]))  # Remove final sugar for 'a'

        fragments["b"].append(prefix_formula)
        fragment_structures["b"].append("--".join(prefix_structure))

        fragments["c"].append(prefix_formula + Counter(structures["inner_link"]))
        fragment_structures["c"].append("--".join(prefix_structure + [format_structure("inner_link")]))

        fragments["d"].append(prefix_formula + Counter(structures["outer_link"]))
        fragment_structures["d"].append("--".join(prefix_structure + [format_structure("outer_link")]))

        # 3' side fragments (w, x, y, z)
        suffix_formula = Counter()
        suffix_structure = []
        for j in range(i, len(sequence)):
            suffix_structure.append(format_structure("sugar"))
            suffix_formula.update(structures["sugar"])
            suffix_structure.append(format_structure(sequence[j]))
            suffix_formula.update(structures[sequence[j]])
            if j > i:  # Add inner/outer linkers for all but the first nucleotide
                suffix_structure.append(format_structure("outer_link"))
                suffix_formula.update(structures["outer_link"])
                suffix_structure.append(format_structure("inner_link"))
                suffix_formula.update(structures["inner_link"])
        fragments["w"].append(suffix_formula - Counter(structures["sugar"]))
        fragment_structures["w"].append("--".join(suffix_structure[:-1]))  # Remove final sugar for 'w'

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
