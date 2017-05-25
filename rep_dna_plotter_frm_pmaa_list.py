__author__ = 'mjohnpayne'
from matplotlib import *
from pylab import *
from numpy import *
from Bio import SeqIO
from Bio import Seq
import subprocess
from time import sleep as sl
import matplotlib.pyplot as plt
import math
import matplotlib.gridspec as gridspec
import scipy


ingenome = SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome.fasta","fasta")

#outpos = open("/Volumes/MP_HD/repDNA_data/Pm_all_genome_rep_counts30mer_1step.txt",'w')


ingene = open("/Users/mjohnpayne/Documents/PhD/pop_paper/pop_ids_rep_loci.txt",'r')
#ingene = open("/Users/mjohnpayne/Documents/PhD/bys/complete_list_bys_ids_f.txt",'r')
#ingene = open("/Volumes/MP_HD/repDNA_data/mannoproteins/gmps_pm_f_expand.txt",'r')


genome = {}

for i in ingenome:
    genome[i.id[:-2]] = i.seq

# concat = ''
#
# for i in ingenome:
#     concat += i.seq
#
#ingenome.close()



# def make_kmer_dict(ingen,stepl,winl):
#     #out_seq_counts = open("/Volumes/MP_HD/repDNA_data/Pm_all_genome_rep_30mer_counts.txt",'w')######### Makes dict/file with each 30 and its count in genome ##############
#     counts = {}
#     for cont in ingen:
#         print cont
#         gen = ingen[cont]
#         for i in range(0,len(gen)-winl,stepl):
#             win = gen[i:i+winl]
#             if win not in counts:
#                 counts[win] = 1
#             else:
#                 counts[win] += 1
#     # for j in counts:
#     #     out_seq_counts.write(str(j) + '\t' + str(counts[j]) + '\n')
#     # out_seq_counts.close()
#     ingenome = SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome.fasta","fasta")######### Makes file with step of 1 and window of 30 accross whole genome with counts of each 30mer in genome ##############
#     counts2 = {}
#     for x in ingenome:
#         print x.id
#         counts2[x.id] = {}
#         gen = str(x.seq)
#         for i in range(0,len(x.seq),stepl):
#             win = gen[i:i+winl]
#             counts2[x.id][str(i)] = counts[win]
#             outpos.write(x.id + '\t' + str(i) + '\t' + str(counts[win]) + '\n')
#     return counts2
#
#
# make_kmer_dict(genome,1,30)

###########################################

# read in file from previous section

insteps = open("/Volumes/MP_HD/repDNA_data/Pm_all_genome_rep_counts30mer_1step.txt",'r')

count_dict = {}

for line in insteps:
    col = line.strip('\n').split('\t')
    if col[0][:-2] not in count_dict:
        print col[0][:-2]
        count_dict[col[0][:-2]] = {col[1]:col[2]}
    else:
        count_dict[col[0][:-2]][col[1]] = col[2]


# outfile.close()
incount = open("/Volumes/MP_HD/repDNA_data/Pm_all_genome_rep_100mer_counts.txt",'r')

count_dict = {}

for line in incount:
    col = line.strip('\n').split('\t')
    count_dict[col[0]] = int(col[1])
print len(count_dict)


def n_random_genes_from_gff(n):
    gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff",'r')
    genes = []
    for line in gff:
        if 'ID=gene' in line:
            col = line.strip('\n').split('\t')
            det = col[8].split(';')
            pmaa = det[1][5:]
            genes.append(pmaa)
    rand_smpl = []
    cnt = 0
    while cnt < n:
        random.seed(cnt+398)### +398 for random set 2
        rand_smpl.append(genes[random.randint(0,len(genes)-1)])
        cnt += 1
    print rand_smpl
    return rand_smpl

def extract_genome_counts_across_gene(gene,genom,prime5,prime3):
    kmer = 100
    gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff",'r')
    seq = ""
    p3 = 0
    p5 = 0
    glen = 0
    for line in gff:
        if 'ID=gene' in line and gene in line:
            col = line.strip('\n').split('\t')
            det = col[8].split(';')
            pmaa = det[1][5:]
            if int(col[3]) > prime5:
                p5 = len(str(genom[col[0]][int(col[3])-prime5:int(col[3])]))
            else:
                p5 = len(str(genom[col[0]][0:int(col[3])]))
            if (int(col[4]) + prime3) > len(genom[col[0]]):
                p3 = len(str(genom[col[0]][int(col[4]):len(genom[col[0]])]))
            else:
                p3 = len(str(genom[col[0]][int(col[4]):int(col[4])+prime3]))
            if int(col[3]) < prime5:
                seq = str(genom[col[0]][0:int(col[4])+prime3])
            elif (int(col[4]) + prime3) > len(genom[col[0]]):
                seq = str(genom[col[0]][int(col[3])-prime5:len(genom[col[0]])])
            else:
                seq = str(genom[col[0]][int(col[3])-prime5:int(col[4])+prime3])
            glen = len(str(genom[col[0]][int(col[3]):int(col[4])]))
    counts = []
    for i in range(len(seq)-kmer):
        win = seq[i:i+kmer]
        if "N" in win or "n" in win:
            counts.append(0)
        else:
            counts.append(count_dict[win])
    gff.close()
    return counts,p5,p3,glen

gcounts = {}
g5primes = {}
g3primes = {}
glens = {}

#ingene = n_random_genes_from_gff(16)

for i in ingene:
    gene = i.strip('\n')
    print gene
    gcounts[gene],g5primes[gene],g3primes[gene],glens[gene] = extract_genome_counts_across_gene(gene,genome,10000,10000)
    ##### converts repetitive counts into 1 no reps and 2 reps #####
    # nlist = []
    # for j in gcounts[gene]:
    #     if j > 1:
    #         nlist.append(2)
    #     else:
    #         nlist.append(1)
    # gcounts[gene] = nlist
    #####
    #print i,gcounts[gene][:50]

cols = 9
rows = 7#int(math.ceil(float(len(gcounts))/cols))
c = 0


fig=plt.figure(figsize=(8.27,11.7))
gs = gridspec.GridSpec(rows,cols,hspace = 0.5)
r = 0
c = 0

for i in list(sorted(gcounts.keys())):
    prime5len = g5primes[i]
    prime3len = g3primes[i]
    genelen = glens[i]
    genpos = [0.1]*prime5len + genelen*[100] + [0.1]*prime3len

    ax1 = plt.subplot(gs[r,c])
    ax0 = plt.subplot(gs[r,c])


    ax1.semilogy(genpos, lw=0.1, antialiased = True, color = 'lightcoral')

    ys = [0.1] + gcounts[i] + [0.1]
    ys = [x+0.1 if x == 0 else x for x in ys]

    ax0.semilogy(ys, lw=0.1, antialiased = True,color = 'steelblue')

    d = scipy.zeros(len(genpos))
    ax1.fill_between([j for j in range(len(genpos))],genpos,where=genpos>=d,interpolate=True,color = 'lightcoral')

    d = scipy.zeros(len(ys))
    ax0.fill_between([j for j in range(len(ys))],ys,where=ys>=d,interpolate=True,color='steelblue')

    ylim(0.8,100)
    xlim(0,len(genpos))
    xticks([len(genpos)/2,len(genpos)],fontsize=8)
    yticks([])
    title(i,fontsize=12)
    if c == 8:
        c = 0
        r += 1
    else:
        c +=1
#
#plt.show()
#
plt.savefig("/Users/mjohnpayne/Documents/PhD/pop_paper/pop_repetitive_loci_smaller.pdf",dpi=200)
# #plt.savefig("/Volumes/MP_HD/repDNA_data/bys_10kb_flanking_repetitiveness_3.pdf",dpi=200)
#plt.savefig("/Volumes/MP_HD/repDNA_data/galmanprot_10kb_flanking_repetitiveness_expand.pdf",dpi=200)
# #plt.savefig("/Volumes/MP_HD/repDNA_data/16random2_10kb_flanking_repetitiveness_2.pdf",dpi=200)