__author__ = 'mjohnpayne'

from matplotlib import *
from pylab import *
from numpy import array
import numpy as np
# from Bio import SeqIO
# from Bio import Seq
# import subprocess
from time import sleep as sl
import matplotlib.pyplot as plt
import math
import matplotlib.gridspec as gridspec
# import scipy

#genlist = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/galactomannoproteins/gmps_pm_f.txt").readlines()
#genlist = open("/Users/mjohnpayne/Documents/PhD/bys/complete_list_bys_ids.txt").read().split('\r')#lines()
#genlist = open("/Users/mjohnpayne/Documents/PhD/ASPS/pop_ids.txt").readlines()
#genlist = open("/Users/mjohnpayne/Documents/PhD/pop_paper/pop_ids.txt").readlines()
#genlist = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/depth_plots_python/control_genes.txt").readlines()
genlist = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/recQ_helicases/recQ_helicase_gene_ids.txt").read().split('\r')#lines()



coverage = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene coverage/cov_comparison_between_strains_corrected_for_tot_reads_for_python.txt",'r').read().split('\r')

strains = ['BR2', 'BR2 SD1', 'BR2 SD2', 'HR2', 'F4', 'FRR3482', 'FRR3840', 'FRR3841', 'FRR3871', 'FRR4059', 'G09027', 'G09027 SD1', 'G09027 SD2', 'G09043', 'G11012', 'G11203', 'G11203 SD1', 'G11203 SD3', 'G11203 SD4', 'G11702']

data = {}
N1 = {}
genes = []

for i in coverage[1:]:
    name = i[0]
    genes.append(name)
    i = i.split('\t')
    data[i[0]] = map(float,i[1:])
    arr = array(map(float,i[1:]))
    arr = arr[arr<160]
    arr = arr[arr>60]
    N1[i[0]] = arr.mean()

strains = coverage[0].strip('\n').split('\t')[1:]
print strains

def n_rand_genes(genelst,no_genes):
    rand_smpl = []
    cnt = 0
    while cnt < n:
        random.seed(cnt+398)### +398 for random set 2
        rand_smpl.append(genes[random.randint(0,len(genes)-1)])
        cnt += 1

maxlst = []
nlist = []

for i in genlist:
    i = i.strip('\n')
    maxlst += [max(data[i])]
    if max(data[i]) > N1[i]*1.5 or min(data[i]) < N1[i]*0.6:
        #print i
        nlist += [i]

print nlist
genlist = nlist

composite_list = [genlist[x:x+40] for x in range(0, len(genlist),40)]

#genlist += ["PMAA_065480"]
#genlist += ["PMAA_007060"]

#print composite_list

maxy = min([600,max(maxlst)])

def plotgen(genlist,number):
    cols = 5
    #rows = int(math.ceil(float(len(genlist))/cols))
    rows = 8

    fig=plt.figure(figsize=(8.27,11.7))
    gs = gridspec.GridSpec(rows,cols,hspace = 0.5)
    r = 0
    c = 0

    for i in genlist:
        print r,c
        ax1 = plt.subplot(gs[r,c])
        ax2 = plt.subplot(gs[r,c])
        ax3 = plt.subplot(gs[r,c])
        ax4 = plt.subplot(gs[r,c])
        ax0 = plt.subplot(gs[r,c])
        covvals = np.array(data[i])
        copy1 = float(N1[i])
        weight = 1.2
        if covvals.max() < copy1*2-40:
            ax1.plot([copy1]*len(strains),lw=weight,color='firebrick')
            ax0.bar(range(len(strains)),covvals,lw=0.1,color = 'steelblue')
        elif covvals.max() < copy1*3-40:
            ax1.plot([copy1]*len(strains),lw=weight,color='firebrick')
            ax2.plot([copy1*2]*len(strains),lw=weight,color='darksage')
            ax0.bar(range(len(strains)),covvals,lw=0.1,color = 'steelblue')
        elif covvals.max() < copy1*4-40:
            ax1.plot([copy1]*len(strains),lw=weight,color='firebrick')
            ax2.plot([copy1*2]*len(strains),lw=weight,color='darksage')
            ax3.plot([copy1*3]*len(strains),lw=weight,color='goldenrod')
            ax0.bar(range(len(strains)),covvals,lw=0.1,color = 'steelblue')
        else:
            ax1.plot([copy1]*len(strains),lw=weight,color='firebrick')
            ax2.plot([copy1*2]*len(strains),lw=weight,color='darksage')
            ax3.plot([copy1*3]*len(strains),lw=weight,color='goldenrod')
            ax4.plot([copy1*4]*len(strains),lw=weight,color='slategrey')
            ax0.bar(range(len(strains)),covvals,lw=0.1,color = 'steelblue')
        #ylim(0,maxy+10)
        xticks([])
        yticks([])

        title(i,fontsize=12)
        if c == cols-1:
            c = 0
            r += 1
        else:
            c +=1

    #plt.show()
    #plt.savefig("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/depth_plots_python/diffGMPs_no_labs.pdf",dpi=300)
    #plt.savefig("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/depth_plots_python/diff_bys_no_labs.pdf",dpi=300)
    #plt.savefig("/Users/mjohnpayne/Documents/PhD/pop_paper/pop_CI_copy_no.pdf",dpi=300)
    #plt.savefig("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/depth_plots_python/control_nolabs.pdf",dpi=300)
    #plt.savefig("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/depth_plots_python/diff_"+str(number)+"_recQ_no_labs.pdf",dpi=300)
    plt.savefig("/Users/mjohnpayne/Documents/PhD/THESIS/Chapter4 - Intraspecies genomics/C4 - figures/"+str(number)+"_recQ_copyno.pdf",dpi=300)


c = 0
for i in composite_list:
    print composite_list[c]
    plotgen(composite_list[c],c)
    c += 1
