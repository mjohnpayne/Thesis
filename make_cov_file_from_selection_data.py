__author__ = 'mjohnpayne'

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import numpy as np
from matplotlib import *
from pylab import *
from numpy import *
from time import sleep as sl

RPKM_index = 15

window = 50000

genome = {}
for record in SeqIO.parse("/Users/mjohnpayne/Documents/Phd/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq

## get gene positions - gff


inpath = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/pos_sel_tests.txt'



def make_cov_dict(selpath):
    pos = {}
    ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')
    rpkms = open(selpath,'r')
    in_rpkms = rpkms.read().split('\n')
    for line in ingff:
        if '#' not in line:
            col = line.strip('\n').split('\t')
            if col[2] == 'gene':
                pmaa = col[8].split(';')[1].replace('Name=','')
                if pmaa[-1] == '0':
                    pos[pmaa] = (float(col[3]),float(col[4]),col[0])

    count = 0
    for line in in_rpkms:
        col = line.strip('\n').split('\t')
        if count == 0:
            #for i in col: print i + ' ' + str(col.index(i))
            count = 1
        else:
            if col[0] in pos:
                pos[col[0]] += (float(col[2]),)

            count +=1

    for i in pos:
        if len(pos[i]) < 4:
            pos[i] += (0,)


    out_dict = {}
    for i in pos:
        if pos[i][2] not in out_dict:
            out_dict[pos[i][2]] = {pos[i][0]:pos[i][3]}
        else:
            out_dict[pos[i][2]][pos[i][0]] = pos[i][3]
        # print pos[i][2],out_dict[pos[i][2]]
        # sl(0.5)
    return out_dict

def make_graph_inp(dict):
    outy = {}
    outx = {}
    for i in dict:
        outy[i] = []
        outpos = sorted(dict[i].keys())
        outx[i] = outpos
        for j in outpos:
            outy[i] += [dict[i][j]]
    return outy,outx


## get RPKM values

##generate output with per position RPKM value


def getavg(val_list):
    avgd = sum(val_list)/len(val_list)
    return avgd

def normalize_val(vals):
    newvals = {}
    for contig in vals:
        newvals[contig] = []
        avg = getavg(vals[contig])
        for x in vals[contig]:
            y = (float(x)/avg)*100
            newvals[contig].append(y)
    return newvals




cov = make_cov_dict(inpath)

yvals,xvals = make_graph_inp(cov)




#position,values = make_cov(cov)
#values = normalize_val(values)

def group_lists(inlist,n):
    outlist = {}
    for j in inlist:
        outlist[j] = []
        for i in range(0,len(inlist[j])-n,n):
            outlist[j] += [np.average(np.asarray(inlist[j][i:i+n]))]
    return outlist

av_over = 1
yvals = group_lists(yvals,av_over)
xvals = group_lists(xvals,av_over)


rel_size = {}


### Finds largest contig

largest = ''
size = 0
for contig in genome:
    if len(genome[contig]) > size:
           size = len(genome[contig])
           largest = contig

##finds relative sizes of all contigs

for contig in genome:
    rel_size[contig] = float((float(len(genome[contig]))/float(len(genome[largest]))))

used_contigs = [x for x in genome.keys() if len(genome[x]) > 50000]

contig_no = (len(used_contigs))
#print rel_size

##makes plot




fig = plt.figure()


counter = 1

for contig in used_contigs:
    avg_depth = getavg(yvals[contig])
    if avg_depth == 0:
        avg_depth = 400
    else:
        avg_depth = avg_depth
    y_min = min(yvals[contig])
    y_max = max(yvals[contig])
    x_start = 0#position[contig][0]
    #pos_len = len(position[contig]) - 1
    x_end = len(genome[contig])#position[contig][pos_len]
    X = numpy.asarray(xvals[contig])#np.linspace(x_start, x_end, pos_len + 1)
    Y = numpy.asarray(yvals[contig])
    new_right = (10 + (85*rel_size[contig]))/100
    gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
    ax0 = plt.subplot(gs[counter-1])
    input1, = ax0.plot(X,Y, lw=0, antialiased = True, alpha=0.5, color = 'black')
    d = np.asarray([0]*len(yvals[contig]))
    #print len(d),len(xvals[contig]),len(yvals[contig])
    ax0.fill_between(xvals[contig], yvals[contig], where=yvals[contig]>=d, interpolate=True, color='red',lw=0.8)
    ax0.fill_between(xvals[contig], yvals[contig], where=yvals[contig]<=d, interpolate=True, color='blue',alpha=1,lw=0)
    # elif len(sys.argv) == 6:
    #     Y2 = numpy.asarray(val2[contig])
    #     input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
    #     Y3 = numpy.asarray(val3[contig])
    #     input3, = ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = name3)
    # elif len(sys.argv) == 7:
    #     Y2 = numpy.asarray(val2[contig])
    #     input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
    #     Y3 = numpy.asarray(val3[contig])
    #     input3, = ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = name3)
    #     Y4 = numpy.asarray(val4[contig])
    #     input4, = ax0.plot(X,Y4, lw=1, antialiased = True, alpha=0.5, color = 'orange', label = name4)
    ylabel(contig, rotation = 'horizontal', fontsize = 10)
    xlim(0,x_end+1)
    ylim(int(y_min),int(y_max))
    xticks([])
    #xticks(range(0,x_end,100000),fontsize = 1)
    yticks([int(y_min),int(y_max)],fontsize = 7)
    counter = counter + 1

# fontP = FontProperties()
# fontP.set_size('small')

#if RPKM_index2 > 0:
#     plt.figlegend([input1,input2], (name1,name2), 'upper right',prop=fontP)
# elif len(sys.argv) == 6:
#     plt.figlegend([input1,input2,input3], (name1,name2,name3), 'upper right',prop=fontP)
# elif len(sys.argv) == 7:
#plt.figlegend([input1,input2,input3,input4], (name1,name2,name3,name4), 'upper right',prop=fontP)
# elif len(sys.argv) == 4:
#plt.figlegend([input1], (name1,), 'upper right',prop=fontP)


#plt.show()
plt.savefig('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/av'+str(av_over)+'_dn-ds_genome_wide_per_gene.pdf',dpi=600,)


