#!/usr/bin/env python

from collections import defaultdict

from math import ceil

## matrix adapted from http://www.compbio.dundee.ac.uk/manuals/amas/amas_manual.txt
matrix = """
ILVCAGMFYWHKREQDNSTP -BZX # AA
1111111111110000001011001 # Hydrophobic
0000000011111111111011111 # Polar
0011110000000001111111001 # Small
0000000000000000000111001 # Proline
0000000000111000000011001 # Positive
0000000000000101000011001 # Negative
0000000000111101000011001 # Charged
0000110000000000010011001 # Tiny
1110000000000000000011001 # Aliphatic
0000000111100000000011001 # Aromatic
""".strip().split("\n")

aa_properties = defaultdict(list)

for i, aa in enumerate(matrix[0][:25]):
    for p in matrix[1:]:
        aa_properties[aa].append( int(p[i]) )


def conservationScore(s):

    aa_count = defaultdict(int)
    n = len(s)

    for aa in s:
        aa_count[aa] += 1

    # if "-" in aa_count:
    #     c = aa_count["-"] - ceil(0.1 * n)

    #     if c <= 0:
    #         del aa_count["-"]
    #     else:
    #         aa_count["-"] = c

    if "-" in aa_count:
        gaps = aa_count["-"]
        del aa_count["-"]
    else:
        gaps = 0

    N = len(aa_count)

    counts = [0]*10

    for aa in aa_count:
        counts = [ (x+y) for (x,y) in zip(counts, aa_properties[aa]) ]

    score = 0

    for count in counts:
        score += (count == 0) or (count == N)

    score = max(score, max(aa_count.values()) * 10 / n)

    return score * (n - gaps) / n

if __name__ == '__main__':
    print conservationScore("DDDDDDDDD")
    print conservationScore("DDDDDDDDE")
    print conservationScore("PPDDDDDE")
    print conservationScore("DDDDEEEEF")
    print conservationScore("DDDDEEEFF")
    print conservationScore("DD-----")
    print conservationScore("DDDDD--")
    print conservationScore("DDDDDD-")


