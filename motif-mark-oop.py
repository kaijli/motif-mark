import os
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


WIDTH, HEIGHT = 1000, 700    # size of canvas
box_height = 25        # height of exon and motif features
start_y = [100, 200, 300, 400, 500, 600]

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

rgba_motifs = [(0.99216,0.90588,0.14510, 1), 
                (0.70980,0.87059,0.16863, 1), 
                (0.43137,0.80784,0.34510, 1), 
                (0.20784,0.71765,0.4745, 1), 
                (0.12157,0.61961,0.53725, 1), 
                (0.14902,0.50980,0.55686, 1), 
                (0.19216,0.40784,0.55686, 1), 
                (0.24314,0.28627,0.53725, 1), 
                (0.28235,0.15686,0.47059, 1), 
                (0.26667,0.00392,0.32941, 1)]



# set display
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
context = cairo.Context(surface)
context.translate(50, 0)

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
    def __init__(self, sequence_, y_pos_):
        '''
        This is the base class for Sequence objects
        '''
        self.sequence = sequence_
        self.length = len(self.sequence)
        self.exon_seq = ""
        self.intron_seq1 = ""
        self.intron2 = ""
        self.y_pos = y_pos_

    def intron_exon(self):
        '''
        This function is applied to Sequence objects to draw the introns and extrons
        found in a gene. 
        '''
        self.exon_seq = re.findall(r'([A-Z]+)', self.sequence)[0]
        exon_start = self.sequence.find(self.exon_seq)
        exon_stop = exon_start + len(self.exon_seq)

        self.intron_seq1 = self.sequence[:exon_start]
        intron_stop1 = exon_start - 1

        self.intron_seq2 = self.sequence[exon_stop:]
        intron_start2 = exon_stop + 1

        self.draw_exon(exon_start, exon_stop)
        self.draw_intron(0, intron_stop1)
        self.draw_intron(intron_start2, self.length)

    def draw_exon(self, exon_start, exon_stop):
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(box_height)
        context.move_to(exon_start, self.y_pos)     # x, y
        context.line_to(exon_stop, self.y_pos)
        context.stroke()
        
    
    def draw_intron(self, intron_start, intron_stop):
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(5)
        context.move_to(intron_start, self.y_pos)     # x, y
        context.line_to(intron_stop, self.y_pos)
        context.stroke()

class Motif(Sequence):
    def __init__(self, sequence_, y_pos_, color_):
        super().__init__(sequence_, y_pos_)
        self.sequence = sequence_
        self.color = color_
        self.start_pos = list()
        self.translated = ""
        self.length = len(self.sequence)
        self.y_pos = y_pos_
 
    def translate(self):
        print(set(self.sequence)<=(ambiguous_bases))
        if set(self.sequence)<=(ambiguous_bases):
            for base in self.sequence:
                self.translated = re.sub(rf'{base}', f"[{bases_represented.get(base)}]", self.sequence)
        else:
            self.translated = self.sequence

    def add_loc(self, sequence: str):
        found_motifs = re.finditer(rf"{self.translated}", sequence.upper())
        self.start_pos = list(x.start() for x in found_motifs)

    def draw(self):
        for start in self.start_pos:
            context.set_source_rgba(self.color[0], self.color[1], self.color[2], self.color[3])
            context.rectangle(start, (self.y_pos-box_height), self.length, (self.y_pos+box_height/2)) 
            context.stroke()
      

oneline_fasta(fasta_file, fasta_oneline)
motif_list = import_motifs(motifs_file)

with open(fasta_oneline, "r") as fr: 
    read_num = 0      
    for line in fr:
        line = line.strip()
        if not line.startswith(">"):
            read = Sequence(line, start_y[read_num])
            read.intron_exon()
            for i in range(len(motif_list)):
                temp_obj = Motif(motif_list[i], start_y[read_num], rgba_motifs[i])
                temp_obj.translate()
                temp_obj.add_loc(read.sequence)
                temp_obj.draw()
            read_num += 1

# draw legend

surface.write_to_png(png_file)

surface.finish()

if os.path.exists(fasta_oneline):
    os.remove(fasta_oneline)

