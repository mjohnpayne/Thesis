import re
from Bio import SeqIO
from Bio import GenBank
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC
import sys
genome ={}
for record in SeqIO.parse('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta', "fasta"):
        genome[record.id] = str(record.seq)

##for record in genome:
##        print record + " " + str(len(genome[record]))
        
insnps = open(sys.argv[1],'r')#open('/Volumes/MP_HD/PhD/clinical_isolate_genomic_data/9_14_PM1_illumina/GATK_processing/PM1_rg_dedup_re_sort.vcf','r')

outfasta = open(sys.argv[2],'a')#open('/Volumes/MP_HD/PhD/clinical_isolate_genomic_data/9_14_PM1_illumina/GATK_processing/PM1_rg_dedup_re_sort.fasta','a')
snps = {}
counter = 0
for line in insnps:
    if '#' not in line:
        columns = line.split("\t")
        if 'CONSVAR' in line and int(columns[7][3:columns[7].find(';AF')]) > 45:
            contig = genome[columns[0]]
            contig = contig[:int(columns[1])] + columns[4] + contig[int(columns[1])+1:]
            genome[columns[0]] = contig
            counter = counter + 1
            if counter%1000 == 0:
                print counter
            # snps[(columns[0],columns[1])] = [columns[0],columns[3],columns[4]]

# for pos in snps:
#     contigid = snps[pos][0]
#     contig = genome[contigid]
#     contig = list(contig)
# #    print contigid
# #    print pos
# #    print contig[int(pos)-1]
#     contig[int(pos[1])-1] = snps[pos][2]
# #    print contig[int(pos)-1]
#     contig = "".join(contig)
#     genome[contigid] = contig
#     counter = counter + 1
#     if counter%100 == 0:
#         print counter

    
for contig in genome:
    contigname = contig
    sequence = str(genome[contig])
    contigname = SeqRecord(Seq(sequence,generic_nucleotide),id=contig)
    SeqIO.write(contigname,outfasta, "fasta")

insnps.close()
outfasta.close()
