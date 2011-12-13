#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import sys
import os
import re

from collections import defaultdict

## The Clustal color scheme, as shown on http://ekhidna.biocenter.helsinki.fi/pfam2/clustal_colours
_rules = """
blue
    (W,L,V,I,M,F): {50%, p}{60%, wlvimafcyhp}
    (A): {50%, p}{60%, wlvimafcyhp}{85%, t,s,g}
    (C): {50%, p}{60%, wlvimafcyhp}{85%, s}
red
    (K,R): {60%, kr}{85%, q}
green
    (T): {50%, ts}{60%, wlvimafcyhp}
    (S): {50%, ts}{80%, wlvimafcyhp}
    (N): {50%, n}{85%, d}
    (Q): {50%, qe}{60%, kr}
pink
    (C): {85%, c}
magenta
    (D): {50%, de,n}
    (E): {50%, de,qe}
orange
    (G): {always}
cyan
    (H,Y): {50%, p}{60%, wlvimafcyhp}
yellow
    (P): {always}
"""

## convert the above rules into a dictionary of functions that
## check how often amino acids occur in a given amino acid set

def get_check(aa_set):
    """
    >>> get_check("a")("abac")
    2
    >>> get_check(None)("abc")
    1
    """
    
    if aa_set is None: return lambda x : 1
    
    r = re.compile(r"[%s]" % aa_set, re.I).findall
    
    return lambda x : len(r(x))

rules = defaultdict(list)
aa_sets = {}

for line in _rules.split("\n"):
    if not line or line.startswith("#"): continue
    if not line.startswith(" "): 
        color = line
        continue
    
    aa, rule1, rule2 = re.search(r"\((.*)\): \{([^}]*)\}(?:\{([^}]*)\})?",line).groups()
    aa = aa.split(",")
    
    for rule in (rule1,rule2):
        if rule is None: continue
        if rule == "always": 
            for a in aa:
                rules[a].append( (0,None,color))
            if None not in aa_sets: aa_sets[None] = get_check(None)
            continue
            
        cutoff = 0.01*int( rule[:rule.index("%")] )
        for aa_set in rule[rule.index(" ")+1:].split(","):
            for a in aa:
                rules[a].append( (cutoff, aa_set, color) )
            if aa_set not in aa_sets: aa_sets[aa_set] = get_check(aa_set)
            
def assign_colors(aa):
    """
    >>> assign_colors("PPP")
    ['yellow', 'yellow', 'yellow']
    >>> assign_colors("PYP")
    ['yellow', 'cyan', 'yellow']
    >>> assign_colors("H")
    ['cyan']
    >>> assign_colors("GHG")
    ['orange', '', 'orange']
    """
    
    checked_sets = {}
    
    colors = []
    
    for a in aa:
        for (cutoff, aa_set, color) in rules[a.upper()]:
            if aa_set not in checked_sets:
                x = checked_sets[aa_set] = aa_sets[aa_set](aa) / len(aa)
            else:
                x = checked_sets[aa_set]

            if x >= cutoff:
                colors.append(color)
                break
            
        else:
            colors.append("")
                
    return colors


if __name__ == '__main__':
    import doctest
    doctest.testmod()



