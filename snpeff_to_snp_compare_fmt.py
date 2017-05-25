__author__ = 'mjohnpayne'

from Bio import SeqIO
from time import sleep as sl
import glob

infile = glob.glob('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/Pass_filter_snps/*_filt.ann.vcf')
in_pmaa = SeqIO.parse('/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/pm_proteins_with_byss.fasta','fasta')


aa = open('/Volumes/MP_HD/SNPeff/Amino_acid_three_to_one.txt','r')
codon_aa = open('/Volumes/MP_HD/SNPeff/codon_to_aa_one_letter.txt','r')

abrev = {}

for line in aa:
    col = line.strip('\n').split('     ')
    abrev[col[0]] = col[1]
codon = {}
for line in codon_aa:
    col = line.strip('\n').split(' ')
    codon[col[0]] = col[1]

pmaas = []
for i in in_pmaa:
    pmaas.append(i.id)

def abbrev_aa(string):
    for ch in abrev.keys():
        if ch in string:
            string=string.replace(ch,abrev[ch])
    return string



strains = []
outdict = {}
for i in infile:
    last = i.split('/')[-1]
    strain = last[:last.find('_')]
    strains.append(strain)
    outfile = open('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/Pass_filter_snps/' + strain + '_filt_SNPeff.txt','w')
    inf = open(i,'r').readlines()
    frameshift_pos = {}
    frameshift_pc_left = {}
    for j in inf:
        if j[0] == '#':
            continue
        elif '|missense_variant|' in j:
            col = j.strip().split('\t')
            inf = col[7].split(';')
            for t in inf:
                if t[:3] == 'ANN':
                    annot = t.split('|')
                    pmaa = annot[3]
                    protchange = abbrev_aa(annot[10])[2:]
                    outfile.write(pmaa + '\t' + protchange + '\n')
    outfile.close()
