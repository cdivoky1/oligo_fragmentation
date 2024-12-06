import argparse

reference = {
    "sugar": {"C":5, "H":9,"O":3},
    "phosphate":{"O":3, "P":1},
    "A":{"C":5, "H":4, "N":5, "O":0},
    "U":{"C":4, "H":4, "N":2, "O":2},
    "C":{"C":4, "H":4, "N":3, "O":1},
    "G":{"C":5, "H":4, "N":3, "O":1},
}

def build_complete(seq):
    # Initialize structure dictionary
    structure = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0}  # Initialize with zeroes for all atoms
    struct_string = []

    for nt in range(0,len(seq)):
        # Add sugar counts
        for atom, count in reference["sugar"].items():
            structure[atom] += count
            
        # Add phosphate counts
        for atom, count in reference["phosphate"].items():
            structure[atom] += count
        # Add base counts
        for atom, count in reference[seq[nt]].items():
            structure[atom] += count
        
    return structure

def generate_a(seq):
   a_frags = []
   frag = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0}  # Initialize with zeroes for all atoms
   
   for nt in range(0,len(seq)):
        for atom, count in reference["sugar"].items():
            frag[atom] += count
        # Add phosphate counts
        for atom, count in reference["phosphate"].items():
            frag[atom] += count
        # Add base counts
        for atom, count in reference[seq[nt]].items():
            frag[atom] += count
        
        print(frag)
        frag["O"] -= 1
        frag["P"] -= 1
        print(frag)
        a_frags.append(frag)
   
    return a_frags[]

def generate_b(formula):
    return formula

def generate_c(formula):
    return formula

def generate_d(formula):
    return formula

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate chemical structure of an RNA sequence.")
    parser.add_argument("sequence", type=str, help="RNA sequence (e.g., AUC).")
    args = parser.parse_args()

    # Extract sequence from command-line arguments
    sequence = args.sequence
    print(f"Sequence to be fragmented: {sequence}")
    print(f"Length of sequence: {len(sequence)}")

    base = build_complete(sequence)
    a_ions = generate_a(base)
    print(a_ions)

if __name__ == "__main__":
    main()