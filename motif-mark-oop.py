import os
import cairo
import argparse
import re

##### User intake #####

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


##### Classes #####

class Sequence:
    '''
    The sequence class describes 1 read of a fasta file or other sequences.
    Its functionalities include deciphering between introns and exons depending on letter capitalization.
    It then draws the introns as lines and exons as boxes. 
    '''
    def __init__(self, sequence_, y_pos_, context_):
        '''
        This is the base class for Sequence objects
        '''
        self.sequence = sequence_
        self.length = len(self.sequence)
        self.exon_seq = ""
        self.intron_seq1 = ""
        self.intron2 = ""
        self.y_pos = y_pos_
        self.context = context_

    def intron_exon(self):
        '''
        This function is applied to Sequence objects to draw the introns and extrons
        found in a gene. 
        There is transparency applied for motifs that may overlap on the reads.
        '''
        self.exon_seq = re.findall(r'([A-Z]+)', self.sequence)[0]
        exon_start = self.sequence.find(self.exon_seq)
        exon_stop = exon_start + len(self.exon_seq)

        self.intron_seq1 = self.sequence[:exon_start]
        intron_stop1 = exon_start

        self.intron_seq2 = self.sequence[exon_stop:]
        intron_start2 = exon_stop

        self.draw_exon(exon_start, exon_stop)
        self.draw_intron(0, intron_stop1)
        self.draw_intron(intron_start2, self.length)

    def draw_exon(self, exon_start, exon_stop):
        '''
        This function draws an exon from start and stop positions relative
        to the beginning of the sequence read. 
        '''
        self.context.set_source_rgba(0, 0, 0, 1)
        self.context.set_line_width(box_height + 10)
        self.context.move_to(exon_start, self.y_pos)     # x, y
        self.context.line_to(exon_stop, self.y_pos)
        self.context.stroke()
        
    
    def draw_intron(self, intron_start, intron_stop):
        '''
        This function draws an intron from start and stop positions relative
        to the beginning of the sequence read. 
        '''
        self.context.set_source_rgba(0, 0, 0, 1)
        self.context.set_line_width(5)
        self.context.move_to(intron_start, self.y_pos)     # x, y
        self.context.line_to(intron_stop, self.y_pos)
        self.context.stroke()

class Motif(Sequence):
    '''
    The motif object is for individual motif sequences to draw unique motifs from user file input
    as different colors on to-scale image of introns and exons. 
    The motifs it can take include those that contain ambiguous bases.
    The motifs with ambiguous bases are translated into regex statements for searching through the main sequence.
    This class inherits features from the sequence class.
    '''
    def __init__(self, sequence_, y_pos_, context_, color_):
        super().__init__(sequence_, y_pos_, context_)
        self.color = color_
        self.start_pos = list()
        self.translated = ""
        self.length = len(self.sequence)
 
    def translate(self):
        '''
        This function translates the sequence of motifs with ambiguous bases
        into a string that contains bracketed options for regex searching later.
        '''
        for base in self.sequence:
            if base in ambiguous_bases:
                self.translated += re.sub(rf'{base}', f"[{bases_represented.get(base)}]", base)
            else:
                self.translated += base

    def add_loc(self, sequence: str):
        '''
        This function locates all the start positions of the motif in the sequence.
        Updates the list of start positions relative to the sequence string.
        '''
        found_motifs = re.finditer(rf"{self.translated}", sequence.upper())
        self.start_pos = list(x.start() for x in found_motifs)

    def draw(self):
        '''
        This function draws motifs as colorful rectangles, with color varying
        between unique motif objects.
        '''
        for start in self.start_pos:
            self.context.set_source_rgba(self.color[0], self.color[1], self.color[2], self.color[3])
            self.context.rectangle(start, (self.y_pos - box_height/2), self.length, box_height)
            self.context.fill()
            self.context.stroke()
      

##### Global functions #####

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


##### Main #####

# process intake files
oneline_fasta(fasta_file, fasta_oneline)
motif_list = import_motifs(motifs_file)

# line count fasta to determine canvas size
line_count = 0
with open(fasta_oneline, 'r') as fp:
    for line_count, line in enumerate(fp):
        pass
    

WIDTH, HEIGHT = 1000, 100*line_count   # size of canvas
box_height = 40        # height of exon and motif features
start_y = range(0, 100*line_count, 100)
# [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

ambiguous_bases = "WSMKRYBDHVNU"
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
                    "N":"ACGTacgt",
                    "U":"Tt"}

rgba_motifs = [(0.90196, 0.62353, 0.00000, .8), # orange
                (0.33725, 0.70588, 0.91373, .8), # sky blue
                (0.00000, 0.61961, 0.45098, .8), # green
                (0.94118, 0.89412, 0.25882, .8), # yellow
                (0.00000, 0.44706, 0.69804, .8), # blue
                (0.83529, 0.36863, 0.00000, .8), # vermillion
                (0.80000, 0.47451, 0.65490, .8)] # pink

# set display
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
context = cairo.Context(surface)
context.set_source_rgba(1, 1, 1, 1)
context.rectangle(0, 0, WIDTH, HEIGHT)  # white background
context.fill()
context.stroke()
context.translate(50, 50)    # shift x by 50 to the right for all following code


# draw legend 
legend_width = 200
legend_height = 11*(len(motif_list)+3)

context.move_to(0, 0)
context.set_source_rgba(0, 0, 0, 1)
context.show_text("Legend")
context.set_line_width(1)
context.rectangle(0, 5, legend_width, legend_height)
context.stroke() 

# add motif labels to legend
for i in range(len(motif_list)):
    # draw colored box
    context.set_source_rgba(rgba_motifs[i][0], rgba_motifs[i][1], rgba_motifs[i][2], 1)
    context.rectangle(5, 10 + (i*11), 10, 10)
    context.fill()
    context.stroke()
    # write motif in black
    context.move_to(20, 20 + (i*11))
    context.set_source_rgba(0, 0, 0, 1)
    context.show_text(f"Motif {motif_list[i]}")

# add exon to legend
context.set_source_rgba(0, 0, 0, 1)
context.rectangle(5, legend_height - 3*(len(motif_list)+3) -2, 10, 10)
context.fill()
context.stroke()
context.move_to(20, 19 + (len(motif_list))*11)
context.show_text("Gene Exon")

# add intron to legend
context.move_to(5, legend_height - (len(motif_list)) -3)
context.set_line_width(3)
context.line_to(15, legend_height - (len(motif_list)) -3)
context.stroke()
context.set_line_width(1)
context.move_to(20, 3 + legend_height - (len(motif_list)+3))
context.show_text("Gene Intron")



context.translate(0, 150)



# parse fasta file and draw 
with open(fasta_oneline, "r") as fr: 
    read_num = 0      
    for line in fr:
        line = line.strip()
        if not line.startswith(">"):
            # draw introns and exons
            read = Sequence(line, start_y[read_num], context)
            read.intron_exon()
            for i in range(len(motif_list)):
                # draw motifs
                temp_obj = Motif(motif_list[i], start_y[read_num], context, rgba_motifs[i])
                temp_obj.translate()
                temp_obj.add_loc(read.sequence)
                temp_obj.draw()
            read_num += 1
        else:
            # label each read with header from FASTA
            context.move_to(0, start_y[read_num] - 2*box_height/3)
            context.set_source_rgba(0, 0, 0, 1)
            context.show_text(f"{line[1:]}")


surface.write_to_png(png_file)
surface.finish()

# remove temp files created
if os.path.exists(fasta_oneline):
    os.remove(fasta_oneline)

