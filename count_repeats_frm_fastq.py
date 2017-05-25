__author__ = 'mjohnpayne'


from Bio import SeqIO
from time import sleep as sl
import glob

seqlist = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/T_flavus/trimmed/*trim*.fastq")#"/Volumes/MP_HD/CI_GENOME_SEQ/CI_trimmed_fastqs/*_trim_*.fastq")
print seqlist
outfile = open("/Volumes/MP_HD/repDNA_data/telomere_reps_counts/Tf_telomere_reps_counts.txt",'w')

def fastq_motif_count(fq):
    inseq = SeqIO.parse(fq,"fastq")
    count = 0
    count2 = 0
    c = 0
    for i in inseq:
        if c%10000 == 0:
            print c
        seq = str(i.seq)
        count2 += seq.count("CCCTAACCCTAA")
        count2 += seq.count("TTAGGGTTAGGG")
        count += seq.count("TCCTAATCCTAA")
        count += seq.count("TTAGGATTAGGA")
        c += 1
    return count, count2

outfile.write("fastq file\tTCCTAATCCTAA count\tCCCTAACCCTAA count\n")

for i in seqlist:
    id = i.split("/")[-1].replace("_trim_","-").strip(".fastq")
    c1,c2 = fastq_motif_count(i)
    outfile.write(id +'\t' + str(c1) + '\t' + str(c2) + '\n')

outfile.close()