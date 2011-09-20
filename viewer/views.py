import django.shortcuts

import sys

import Bio
from Bio import SeqIO

from collections import defaultdict

import amview.viewer.alignment as alignment

import newick
import Taxonomy

from sqlalchemy import create_engine

tax = Taxonomy.Taxonomy(create_engine('sqlite:///ncbi_taxonomy.db'), Taxonomy.ncbi.ranks)

class CollectLeafs(newick.tree.TreeVisitor):
    
    def __init__(self):
        self.leafs = [[]]
    
    def pre_visit_tree(self,t):
        self.leafs.append([])

    def post_visit_tree(self,t):
        l = self.leafs.pop()
        if "9606" in l:
            l.extend(self.leafs[-1])
            self.leafs[-1] = l
        else:
            self.leafs[-1].extend( l )
    
    def visit_leaf(self, l):
        self.leafs[-1].append(l.identifier)
        
    def get_leafs(self):
        return self.leafs[0]
    

def index(request):
    
    first_col = []
    species_col = []
    body = []

    sequences = []
    rows = []    
    
    _anns = {}
    # for record in SeqIO.parse(open("examples/kif2c_coils.fasta"), "fasta"):
    #     _anns[record.id] = record.seq

    # for line in open("/Volumes/cellapps.biotec.tu-dresden.de/data/cellnet/michaelk/data/ccaligner/projects/T07C4.10/full.pc.tsv"): 
    # for line in open("/Volumes/cellapps.biotec.tu-dresden.de/data/cellnet/michaelk/data/ccaligner/projects/sas-4/pruned.pc.tsv"): 
    for line in open("/Volumes/cellapps.biotec.tu-dresden.de/data/cellnet/michaelk/data/ccaligner/projects/spd-5/amview/full.pc.tsv"): 
    # for line in open("examples/pcm1.tsv"):
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            l = []
            _anns[ line[1:] ] = l
        else:
            (_, register, pvalue) = line.split("\t")
            if float(pvalue) > 0.1: register = "-"
            l.append(register)

    anns = []

    items_per_species = defaultdict(list)

    # for record in SeqIO.parse(open("/Users/mkuhn/Desktop/spd-5/spd5_aligned.faa"), "fasta"):
    # for record in SeqIO.parse(open("/Volumes/cellapps.biotec.tu-dresden.de/data/cellnet/michaelk/data/ccaligner/projects/T07C4.10/full.faa.best"), "fasta"):
    for record in SeqIO.parse(open("/Volumes/cellapps.biotec.tu-dresden.de/data/cellnet/michaelk/data/ccaligner/projects/spd-5/amview/full.faa.best"), "fasta"):
    # for record in SeqIO.parse(open("/Volumes/cellapps.biotec.tu-dresden.de/data/cellnet/michaelk/data/ccaligner/projects/sas-4/pruned.faa.best"), "fasta"):
    # for record in SeqIO.parse(open("/Volumes/cellapps.biotec.tu-dresden.de/data/cellnet/michaelk/data/ccaligner/projects/CEP290/full.faa.best"), "fasta"):
    # for record in SeqIO.parse(open("examples/pcm1_metazoa.faa"), "fasta"):
        species = record.id.split(".", 1)[0]
        items_per_species[species].append(record)
        

    species_tree = tax.tree_lineage(items_per_species.keys())
    
    # go through tree and collect leaf identifiers
    
    tv = CollectLeafs()
    species_tree.dfs_traverse(tv)

    for species in tv.get_leafs():
        
        for record in items_per_species[species]:
            first_col.append("""%s<br />""" % record.id)
            _species = tax.primary_from_id(species)
            species_col.append("""<a href="http://en.wikipedia.org/wiki/Special:Search?search=%s" target="_blank">%s</a><br />""" % (_species,_species))
            anns.append(_anns[record.id])
            sequences.append(record.seq)
            rows.append([])
        
    ann_position = [0]*len(sequences)
        
    for x in range(len(sequences[0])):
        aa = "".join(sequence[x] for sequence in sequences)
        colors = alignment.assign_colors(aa)
            
        for y,(row,a,c) in enumerate(zip(rows,aa,colors)):
            
            classes = []
            
            if c: classes.append(c[0])
            
            if a != "-":
                ann = anns[y][ann_position[y]]
                if ann not in (".-"): classes.append("a"+ann)
                ann_position[y] += 1
            
            if classes:
                row.append("<span class='%s'>%s</span>" % (" ".join(classes),a))
            else:
                row.append(a)
                # row.append("<td>%s</td>" % (a,))

    for row in rows:
        body.append("".join(row)+"<br/>\n")
        # body.append("<tr>"+"".join(row)+"</tr>\n")
    
    ncols = len(rows[0])
    
    table_first_row = "\n".join( ['<div class="hd1">1<br/>|</div>'] + [ '<div class="hd2">%d<br/>|</div>' % i for i in range(5, ncols, 5) ])
    

    table_first_col = "\n".join(first_col)
    table_body = "".join(body)
    
    species_col = "\n".join(species_col)
    
    return django.shortcuts.render_to_response('viewer.html',locals())


