# oligo_fragmentation
Fragments a given RNA (A/C/G/U) oligo into a, b, c, d, w, x, y, and z fragments in miRNA architecture (terminal 5' phosphate and 3' OH). 

Fragmentation breakdown:
![image](https://github.com/user-attachments/assets/ed8481eb-fc1c-4689-be63-43a91642a46f)

Image ref: https://link.springer.com/article/10.1140/epjd/e2011-20616-y (fig 22)

Notes:
- While it prints the reverse ion fragments (w,x,y, & z) in decending order (fragment N -> fragment 1), the ion fragments are stored in their data structure in ascending order (fragment 1 -> fragment N).

Upcoming Versions:
- Handle fragmentation of digested RNA species
- Handle non-miRNA structured RNAs
- Peptide fragmentation 
