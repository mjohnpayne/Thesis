__author__ = 'mjohnpayne'

import matplotlib.pyplot as plt


def filt_in(invcf):
    outls = []
    sizes = []
    pesupport = []
    srsupport = []
    inf = open(invcf,"r")
    for i in inf:
        if i[0] != "#":
            col = i.split("\t")
            details = col[7].split(';')
            pe_reads = int(details[8].replace('PE=',''))
            pesupport +=[pe_reads]
            st = int(col[1])
            end = int(details[6][4:])
            endcont = details[5][5:]
            size = end - st
            sizes += [size]
            if details[0] == "PRECISE":
                sup_reads = int(details[10].replace('SR=',''))
                srsupport += [sup_reads]
    return sizes,pesupport,srsupport

s,p,sr = filt_in("/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filt_precise_all_DEL.vcf")

plt.hist(s, bins=100,color="teal")
plt.xlim([0,max(s)])
plt.show()

plt.hist(p, bins=100,color="teal")
plt.xlim([0,max(p)])
plt.show()

plt.hist(sr, bins=100,color="teal")
plt.xlim([0,max(sr)])
plt.show()