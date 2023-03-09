import os
import numpy as np
import cairo
import argparse
import re

# lookahead function?
##### Global variables and functions #####

def get_args():
    parser = argparse.ArgumentParser(description="A program to draw sequence diagram with introns, exons, and motifs.")
    parser.add_argument("-f", "--fasta", help="fasta file name")
    parser.add_argument("-m", "--motifs", help="motifs file name")
    return parser.parse_args()

args = get_args()
fasta_file = args.fasta
motifs_file = args.motifs
png_file = re.sub(r'(.fasta$)', ".png", fasta_file)
fasta_oneline = re.sub(r'(.fasta$)', "_oneline.fasta", fasta_file)


WIDTH, HEIGHT = 500, 500    # size of canvas


DNA_bases = set("ACGTNacgtn") #each individual base as capital or lowercase
RNAbases = set('ACGUNacgun')
ambiguous_bases = set("WSMKRYBDHVN")
bases_represented = {"W":"ATat",
                    "S":"CGcg",
                    "M":"ACac",
                    "K":"GTgt",
                    "R":"AGag",
                    "Y":"CTct",
                    "B":"CGTcgt",
                    "D":"AGTagt",
                    "H":"ACTact",
                    "V":"ACGacg",
                    "N":"ACGTacgt"}

rgba_motifs = [(253, 231, 37, 1), 
                (181, 222, 43, 1), 
                (110, 206, 88, 1), 
                (53, 183, 121, 1), 
                (31, 158, 137, 1), 
                (38, 130, 142, 1), 
                (49, 104, 142, 1), 
                (62, 73, 137, 1), 
                (72, 40, 120, 1), 
                (68, 1, 84, 1)]

box_height = 100        # height of exon and motif features

# set display
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
context = cairo.Context(surface)


def oneline_fasta(filein: str, fileout: str):
    """
    This function takes a fasta file where the sequences have line breaks in them.
    It takes two filenames, and makes a new files with the name input second. 
    """
    collect_seqs = {}
    with open(filein, "r") as fh:
        header = ""
        for line in fh:
            line = line.strip()
            if ">" in line:
                header = line
            else:
                if header not in collect_seqs.keys():
                    collect_seqs[header] = line
                else:
                    collect_seqs[header] += line
    with open(fileout, "w") as fh:
        for header, seq in collect_seqs.items():
            fh.write(f"{header}\n{seq}\n")

def import_motifs(motif_file: str):
    """
    This function takes a motif file with one motif per line.
    It returns a list of motifs.
    Not case sensitive
    """
    collect_motifs = list()
    with open(motif_file, "r") as fh:
        motif_counter = 0
        for line in fh:
            line = line.strip()
            collect_motifs.append(line.upper())
    return collect_motifs


##### Classes #####

class Sequence:
    def __init__(self, sequence_):
        '''
        This is the base class for Sequence objects
        '''
        self.sequence = sequence_
        self.length = len(self.sequence)

    def intron_exon(self):
        '''
        This function is applied to Sequence objects to draw the introns and extrons
        found in a gene. 
        '''
        exon_seq = re.findall(r'([A-Z]+)', self.sequence)[0]
        exon_start = self.sequence.find(exon_seq)
        exon_stop = exon_start + len(exon_seq)
        intron_seq1 = self.sequence[:exon_start]
        intron_stop1 = exon_start - 1
        intron_seq2 = self.sequence[exon_stop:]
        intron_start2 = exon_stop + 1

        self.draw_exon(exon_start, exon_stop)
        self.draw_intron(0, intron_stop1)
        self.draw_intron(intron_start2, self.length)

    def draw_exon(self, start_, stop_):

        pass
    
    def draw_intron(self, start_, stop_):
        pass

class Motif():
    def __init__(self, sequence_, color_):
        self.sequence = sequence_
        self.color = color_
        self.start_pos = list()
        self.translated = ""

    def translate(self):
        for base in self.sequence:
            if base in ambiguous_bases:
                self.translated = re.sub(rf'{base}', f"[{bases_represented.get(base)}]", self.sequence)
        return self.translated
    
    def add_loc(self, sequence: str):
        found_motifs = re.finditer(rf"{self.translated}", sequence.upper())
        self.start_pos = [x.start() for x in found_motifs]
    def draw(self):
        pass


oneline_fasta(fasta_file, fasta_oneline)
motif_list = import_motifs(motifs_file)

with open(fasta_oneline, "r") as fr:       
    for line in fr:
        line = line.strip()
        if not line.startswith(">"):
            # collect gene name
            pass
        else:
            read = Sequence(0, len(line), line)
            read.intron_exon()
            for motif in motif_list:
                found_motifs = re.finditer(rf"{motif}", line.upper())
                [x.start() for x in found_motifs]


if os.path.exists(fasta_oneline):
    os.remove(fasta_oneline)
    print(f"removed {fasta_oneline}")

