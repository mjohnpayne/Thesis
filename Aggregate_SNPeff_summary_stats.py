__author__ = 'mjohnpayne'

import glob
from time import sleep as sl
import html2text


instats = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_snps/*.html")

# outfile = open("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_snps/snp_stats.txt","w")
#
# outfile.write("strain\tno_snps\tno_missense\tno_nonsense\tno_silent\n")

for i in instats:
    spec = i.split('/')[-1].replace("_GATK_filtered_snps_pass_stats.html","")
#    outfile.write(spec+"\t")
    inf = open(i,"r").read()
    f = html2text.html2text(inf).split('\n')
    #f = open(i,'r')
    print i,spec
    # for j in f[:200]:
    #     print j
    # sl(15)
    for j in f:
        if "** SNP **" in j:
            no = j.replace("** SNP ** |  ","").strip('\n')
#            outfile.write(no+"\t")
        elif "** MISSENSE **" in j:
            no = j.replace("** MISSENSE ** |    |  ","")
            no = no[:no.find("|")-2]
#            outfile.write(no+"\t")
        elif "** NONSENSE **" in j:
            no = j.replace("** NONSENSE ** |    |  ","")
            no = no[:no.find("|")-2]
#            outfile.write(no+"\t")
        elif "** SILENT ** " in j:
            no = j.replace("** SILENT ** |    |  ","")
            no = no[:no.find("|")-2]
#            outfile.write(no+"\n")
#outfile.close()