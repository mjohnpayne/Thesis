# -*- coding: utf-8 -*-
import random
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC
from matplotlib import *
from pylab import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties

insnps = '/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_snps/012_GATK_filtered_snps_pass.ann.vcf'#sys.argv[1]



genome = {}
for record in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
       genome[record.id] = record.seq
       print record.id + "\t" + str(len(genome[record.id]))
#print genome



def make_windows_frm_vcf(insnp):
    #### process vcf gives dict of contig = [positions]
    snps = {}
    contig = ""
    vcf_open = open(insnp,"r")
    for line in vcf_open:
        if '#' not in line:
            col = line.split('\t')
            if contig == col[0]:
                snps[col[0]].append(int(col[1]))
            else:
                contig = col[0]
                snps[col[0]] = [int(col[1])]

    ### generate list of tuples with contig, winstart, winend, snp_count
    winstart = 0
    winend = 10000
    step = 1000
    windowcount = 0
    index = 0
    outls = []

    for contig in snps:
       winstart = 0
       winend = 10000
       step = 1000
       if len(genome[contig]) < winend:
           count = len(snps[contig])
       else:
           while winend < len(genome[contig]):
               while index < len(snps[contig]):
                   if snps[contig][index] > winstart and snps[contig][index] < winend:
                       windowcount = windowcount + 1
                       index = index + 1
                   else:
                       index = index + 1
               outls.append((contig,str(winstart),str(winend),str(windowcount)))
               windowcount = 0
               winstart = winstart + step
               winend = winend + step
               index = 0


    vcf_open.close()





    ######



    ## generate corresponding dictionaries of window positions and window values

    contigs = []
    values = {}
    position = {}

    for columns in outls:
        if columns[0] in contigs:
            values[columns[0]] = values[columns[0]] + [columns[3]]
            position[columns[0]] =  position[columns[0]] + [(int(columns[1])+5000)/1000]
        else:
            contigs = contigs + [columns[0]]
            values[columns[0]] = [columns[3]]
            position[columns[0]] = [(int(columns[1])+5000)/1000]

    newpos = {}
    newval = {}

    for contig in position:
        x_start = position[contig][0]
        pos_len = len(position[contig]) - 1
        x_end = position[contig][pos_len]
        if x_end > 10:
            newpos[contig] = position[contig]
            newval[contig] = values[contig]

    position = newpos
    values = newval
    return contigs,values,position





contigs,values,position = make_windows_frm_vcf(insnps)

rel_size = {}

largest = ''
size = 0
for contig in position:
    if len(position[contig]) > size:
           size = len(position[contig])
           largest = contig




for contig in position:
    rel_size[contig] = float((float(len(position[contig]))/float(len(position[largest]))))

contig_no = (len(position))



fig = plt.figure


counter = 1



for contig in position:
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    X = np.linspace(x_start, x_end, pos_len + 1)
    Y = numpy.asarray(values[contig])
    new_right = (10 + (85*rel_size[contig]))/100
    gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
    ax0 = plt.subplot(gs[counter-1])
    ax0.plot(X,Y, lw=1.0, antialiased = True, alpha=0.9, color = 'red')
    ylabel(contig, rotation = 'horizontal', fontsize = 12)
    xlim(0,x_end+1)
    ylim(0,50)
    xticks([x_end],fontsize = 7), yticks(fontsize = 6)
    counter = counter + 1


plt.show()
# outpdf = sys.argv[2]
# plt.savefig(outpdf,dpi=800)



            
