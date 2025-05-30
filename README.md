# Fragmention Scripts

Fragments a given sequence into a/a-nh/b/c/d/w/x/y/z fragments for oligos and a/b/c/x/y/z fragments for peptides. 

---

### Oligo Fragmentation Pattern

![Oligo Fragmentation Schematic](https://github.com/user-attachments/assets/ed8481eb-fc1c-4689-be63-43a91642a46f)

*Image adopted from [SpringerLink, Eur. Phys. J. D (Figure 22)](https://link.springer.com/article/10.1140/epjd/e2011-20616-y)*

---

### Peptide Fragmentation Pattern

![Peptide Fragmentation Schematic](https://www.matrixscience.com/images/cleavages.gif)

*Image adopted from [Matrix Science](https://www.matrixscience.com/help/fragmentation_help.html)*

---

Notes:
- While it prints the reverse ion fragments in descending order (fragment N → fragment 1), the ion fragments are stored in their data structure in ascending order (fragment 1 → fragment N).
- prints to the screen and doesn't write to a file
