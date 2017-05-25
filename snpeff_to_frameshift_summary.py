__author__ = 'mjohnpayne'

from Bio import SeqIO
from time import sleep as sl
import glob

infile = glob.glob('/Volumes/MP_HD/bys_snpeff_catchup/initial_snpeff_indel_vcfs/*_bys_catchup.vcf')
print infile
in_pmaa = SeqIO.parse('/Volumes/MP_HD/bys_snpeff_catchup/catchup_bys_prots.fasta','fasta')

outfile = open('/Volumes/MP_HD/bys_snpeff_catchup/frameshifts_bys_catchup.txt','w')

pmaas = []
for i in in_pmaa:
    pmaas.append(i.id)

strains = []
outdict = {}
for i in infile:
    last = i.split('/')[-1]
    strain = last[:last.find('_')]
    strains.append(strain)
    inf = open(i,'r').readlines()
    frameshift_pos = {}
    frameshift_pc_left = {}
    for j in inf:
        if j[0] == '#':
            continue
        elif 'frameshift_variant' in j:
            col = j.strip().split('\t')
            inf = col[7].split(';')
            for t in inf:
                if t[:3] == 'ANN':
                    annot = t.split('|')
                    pmaa = annot[3]
                    pos = int(annot[13][:annot[13].find('/')])
                    size = int(annot[13][annot[13].find('/')+1:])
                    perc = 100-(float(pos)/size*100)
                    frameshift_pos[pmaa] = pos
                    frameshift_pc_left[pmaa] = perc
    for k in pmaas:
        if k not in outdict:
            outdict[k] = {}
            if k in frameshift_pos:
                outdict[k][strain] = frameshift_pc_left[k]
            else:
                outdict[k][strain] = '-'
        else:
            if k in frameshift_pos:
                outdict[k][strain] = frameshift_pc_left[k]
            else:
                outdict[k][strain] = '-'


outfile.write('GeneID' + '\t' + '\t'.join(sorted(strains)) + '\n')


for j in sorted(pmaas):
    outfile.write(j + '\t')
    for k in sorted(strains):
        outfile.write(str(outdict[j][k]) + '\t')
    outfile.write('\n')

outfile.close()