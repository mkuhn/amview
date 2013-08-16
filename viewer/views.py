from __future__ import division

import django.shortcuts

import os
import sys
import re

from collections import defaultdict

import amview.viewer.alignment as alignment
import settings

import newick
import Taxonomy

import urllib

from conservation import conservationScore

from sqlalchemy import create_engine

assert os.path.exists("ncbi_taxonomy.db") and os.path.getsize("ncbi_taxonomy.db") > 0, "Could not find NCBI Taxonomy database"

tax = Taxonomy.Taxonomy(create_engine('sqlite:///ncbi_taxonomy.db'), Taxonomy.ncbi.ranks)

class Record(object):
    """Thin clone of BioPython's Sequence record"""
    def __init__(self, identifier, seq):
        self.id = identifier
        self.seq = seq

def parseFasta(fh):
    """A function to extract sequence records from a FASTA file to avoid a BioPython dependency."""

    record_seq = []
    record_id = None

    for line in fh:
        line = line.strip("\n")

        if line.startswith(">"):

            if record_seq:
                yield Record(record_id, "".join(record_seq))

            record_id = line[1:].split()[0]
            record_seq = []
        else:
            record_seq.append(line)

    if record_seq:
        yield Record(record_id, "".join(record_seq))


class CollectLeafs(newick.tree.TreeVisitor):
    """Newick tree visitor that projects a tree onto a list with human at the top."""

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


def index(request, path):

    assert path.endswith("/")

    min_rows = int(request.REQUEST.get("min_rows", 0))
    svg = int(request.REQUEST.get("svg", 0))

    first_col = []
    species_col = []
    body = []

    sequences = []
    rows = []

    _residue_annotations = defaultdict(list)

    # try to load all annotation files
    for filename in settings.ANNOTATION_FILES:

        filename = settings.ALIGNMENT_PATH + path + filename
        if not os.path.exists(filename): continue

        if "gps.txt" in filename:
            load_gps_annotation_file(_residue_annotations, filename)
        else:
            load_standard_annotation_file(_residue_annotations, filename)


    items_per_species = defaultdict(list)

    # load only one multiple alignment
    for filename in settings.ALIGNMENT_FILES:

        filename = settings.ALIGNMENT_PATH + path + filename
        if not os.path.exists(filename): continue

        for record in parseFasta(open(filename)):
            if re.match(r"\d+\.", record.id):
                species = record.id.split(".", 1)[0]
                items_per_species[species].append(record)
            else:
                items_per_species[""].append(record)

        break

    else:
        print >> sys.stderr, "Not found", filename


    assert items_per_species
    species_tree = tax.tree_lineage( k for k in items_per_species.keys() if k )

    # go through tree and collect leaf identifiers

    species_list = []

    if not all( s == "" for s in items_per_species.keys() ):
        tv = CollectLeafs()
        species_tree.dfs_traverse(tv)
        species_list.extend(tv.get_leafs())

    species_list.append("")

    residue_annotations = []
    record_ids = []

    for species in species_list:

        for record in sorted( items_per_species[species], key=lambda(x) : len(x.seq) ):

            record_ids.append(record.id)
            first_col.append("""<span class='record_id'>%s</span><br />""" % record.id)
            if species:
                _species = tax.primary_from_id(species)
                __species = urllib.quote(_species)
                species_col.append("""<a href="http://en.wikipedia.org/wiki/Special:Search?search=%s" target="_blank">%s</a><br />""" % (__species,_species))
            if record.id in _residue_annotations:
                residue_annotations.append(_residue_annotations[record.id])
            else:
                print >> sys.stderr, "No annotation for", record.id

                residue_annotations.append(None)
            sequences.append(record.seq)
            rows.append([])

    ann_position = [0]*len(sequences)

    domain_annotations = {}
    domain_annotations["6239.F56A3.4"] = (
        ("air-1", "lightblue", 805//3, 1410//3),
        ("spd-2", "lightgreen", 604//3, 1008//3),
        ("rsa-2", "lightcoral", 1207//3, 1410//3),
        ("S", "red", 25 ,25 ),
        ("S", "red", 201,201),
        ("S", "red", 387,387),
        ("S", "red", 530,530),
        ("S", "red", 635,635),
        ("S", "red", 649,649),
        ("S", "red", 658,658),
        ("S", "red", 667,667),
        ("S", "red", 697,697),
        ("S", "red", 736,736),
    )

    overview = []

    if svg:
        ncols = render_svg_alignment(body, rows, domain_annotations, min_rows, record_ids, residue_annotations, sequences, ann_position)

        l = []
        l.append("""<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width=800 height=20 preserveAspectRatio="none">""")

        step = 200
        for i in range(step, ncols, step):
            x = int(i/ncols*800.0)
            l.append('<text x="%d" y=15>%d</text>' % (x, i))
            l.append("""<line x1="%d" y1="15" x2="%d" y2="20" style="stroke:rgb(0,0,0);stroke-width:2"/>""" % (x,x))

        l.append("""</svg>""")

        table_first_row = "\n".join(l)

    else:
        cons_scores = render_text_alignment(body, rows, domain_annotations, min_rows, record_ids, residue_annotations, sequences, ann_position)
        ncols = len(rows[0])
        table_first_row = "\n".join( ['<div class="hd1">1<br/>|</div>'] + [ '<div class="hd2">%d<br/>|</div>' % i for i in range(5, ncols, 5) ])
        render_svg_conservation(overview, cons_scores)
        ncol = len(cons_scores)


    table_first_col = "\n".join(first_col)
    table_body = "".join(body)
    species_col = "\n".join(species_col)
    overview = "\n".join(overview)

    return django.shortcuts.render_to_response('viewer.html', locals())

def load_gps_annotation_file(_residue_annotations, filename):

    fh_in = open(filename)
    fh_in.next()

    for line in fh_in:
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            l = defaultdict(lambda : "-")
            _residue_annotations[ line[1:] ].append(l)
        else:
            # 36      LDVISDTSGLGNGVL 2.681   2.59    Phosphorylation
            (residue, _, score, _, _) = line.split("\t")
            if float(score) > 4:
                l[int(residue)-1] = "P"

def load_standard_annotation_file(_residue_annotations, filename):

    for line in open(filename):
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            l = []
            _residue_annotations[ line[1:] ].append(l)
        else:
            if "\t" in line:
                (_, register, pvalue) = line.split("\t")
                if float(pvalue) > 0.1: register = "-"
                l.append(register)
            else:
                l.extend(line)

def getSVG(element, a, b):

    s = "stroke:rgb(0,0,0);stroke-width:2"

    if element == "CC":
        s = "stroke:rgb(0,185,0);stroke-width:6"

    elif element == "P":
        assert a == b
        return """<path d ="M %d 6 L %d 1 L %d 1" fill="rgb(255,0,255)"/>""" % (a,a-5,a+5)

    return """<line x1="%d" y1="5" x2="%d" y2="5" style="%s"/>""" % (a,b,s)


def render_svg_alignment(body, rows, domain_annotations, min_rows, record_ids, residue_annotations, sequences, ann_position):

    element_per_row = {}

    draw_on_top_rows = [ [] for _ in rows ]

    cols = 0

    for x in range(max(len(sequence) for sequence in sequences)):

        aa = "".join( (sequence[x] if x < len(sequence) else "-") for sequence in sequences)

        if min_rows > 0 and len(aa) - aa.count("-") < min_rows: continue

        cols += 1

        for y,(row,draw_on_top,a) in enumerate(zip(rows,draw_on_top_rows,aa)):

            current_ann_position = ann_position[y]

            last_element, last_start = element_per_row.get(y, (None, None))

            current_element = None

            if a != "-":
                current_element = "x"
                for annotations in residue_annotations[y]:
                    ann = annotations[current_ann_position]
                    if ann is not None:
                        if ann in "abcdefg":
                            current_element = "CC"
                        if ann == "P":
                            draw_on_top.append( getSVG("P", x-1, x-1))

                ann_position[y] += 1

            if current_element != last_element:
                if last_element is not None:
                    row.append( getSVG(last_element, last_start, x-1))
                element_per_row[y] = (current_element, x)



    for row, draw_on_top in zip(rows, draw_on_top_rows):
        body.append("""<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width=800 height=10 viewBox="0 0 %d 10" preserveAspectRatio="none">""" % x)
        body.append("".join(row))
        body.append("".join(draw_on_top))
        body.append("""</svg><br/>\n""")

    return cols



def render_svg_conservation(body, cons_scores):

    n = len(cons_scores)
    body.append("""<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width=100%% height=20px viewBox="0 0 %d 20" preserveAspectRatio="none" id="overviewcanvas">""" % n)
    body.append("""  <defs>
        <style type="text/css"><![CDATA[
          line {
            stroke: #333;
            stroke-width: 1
          }
        ]]></style>
  </defs>""")


    body.append("""<line x1="0" y1="0" x2="%d" y2="0" style="stroke: #777"/>""" % n)
    body.append("""<line x1="0" y1="10" x2="%d" y2="10" style="stroke: #ccc"/>""" % n)

    for x, c in enumerate(cons_scores):
        if c:
            body.append("""<line x1="%d" y1="%d" x2="%d" y2="20"/>""" %(x,20-2*c,x))

    body.append("""</svg>\n""")


def render_text_alignment(body, rows, domain_annotations, min_rows, record_ids, residue_annotations, sequences, ann_position):
    annotation_rows = defaultdict(list)
    active_rows = defaultdict(dict)
    cons_scores = []

    for x in range(len(sequences[0])):

        aa = "".join(sequence[x] for sequence in sequences)

        skip_this_column = min_rows > 0 and len(aa) - aa.count("-") < min_rows

        if not skip_this_column:
            cons_scores.append( conservationScore(aa) )

        colors = alignment.assign_colors(aa)

        for y,(row,a,c) in enumerate(zip(rows,aa,colors)):

            record_id = record_ids[y]
            current_ann_position = ann_position[y]

            for domain_id, (domain_name, domain_color, domain_start, domain_end) in enumerate(domain_annotations.get( record_id, () )):

                domain_padding_end = domain_end + 5

                if domain_start <= current_ann_position <= domain_padding_end:

                    if domain_id not in active_rows[record_id]:
                        ## find a new inactive row
                        for current_row in range(0, len(active_rows[record_id])+1):
                            if current_row not in active_rows[record_id].values(): break
                        active_rows[record_id][domain_id] = current_row
                    else:
                        current_row = active_rows[record_id][domain_id]

                    # NOTE: this will truncate the name if it's longer than the domain (plus padding)
                    # NOTE: this ignores gaps

                    offset = current_ann_position - domain_start
                    cell = []

                    if a != "-":
                        if offset == 0:
                            # add new row, if needed
                            if current_row >= len(annotation_rows[record_id]):
                                assert current_row == len(annotation_rows[record_id])
                                annotation_rows[record_id].append([])

                            # add padding to get the starting position right
                            if len(annotation_rows[record_id][current_row]) < (x-1):
                                annotation_rows[record_id][current_row].extend( ["&nbsp;"] * (x - len(annotation_rows[record_id][current_row]) -1) )

                            cell.append("<span style='background-color: %s'>" % domain_color)

                        # print name or fill with space
                        if offset < len(domain_name):
                            cell.append(domain_name[offset])
                        else:
                            cell.append("&nbsp;")

                        # close the span at the domain end
                        if current_ann_position == domain_end:
                            cell.append("</span>")

                        # give back the active row
                        if current_ann_position == domain_padding_end:
                            active_rows[record_id].pop(domain_id)

                    elif offset > 0:
                        cell.append("&nbsp;")

                    if cell:
                        annotation_rows[record_id][current_row].append("".join(cell))

            classes = []
            if c: classes.append(c[0])

            if a != "-":

                if residue_annotations[y] != None:
                    for annotations in residue_annotations[y]:
                        ann = annotations[current_ann_position]
                        if ann not in (".-"): classes.append("a"+ann)

                ann_position[y] += 1

            if skip_this_column:
                continue

            if classes:
                row.append("<span class='%s'>%s</span>" % (" ".join(classes),a))
            else:
                row.append(a)
                # row.append("<td>%s</td>" % (a,))

    for record_id in record_ids:
        for annotation_row in annotation_rows.get(record_id, ()):
            rows.append(annotation_row)

    row = []
    for c in cons_scores:
        a = "-"
#        row.append("<span class='c%d'>%s</span>" % (c,a))
        row.append("<span>%s</span>" % (c if c < 10 else ":"))

    rows.append(row)

    for row in rows:
        body.append("".join(row)+"<br/>\n")
        # body.append("<tr>"+"".join(row)+"</tr>\n")

    return cons_scores

