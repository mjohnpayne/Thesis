__author__ = 'mjohnpayne'
import subprocess
from Bio import SeqIO
import os
import glob
from time import sleep as sl

# infiles = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/*.fasta")
#
# done = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/align_done.txt"
# for i in infiles:
#     d = open(done,"r").readlines()
#     d = [x.strip('\n') for x in d]
#     print float(len(d)*100/len(infiles))
#     pmaa = i.split('/')[-1][:11]
#     if pmaa not in d:
#         mega_align_args = "/usr/local/bin/megacc -a /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/clustal_align_coding.mao -o /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/alignments -s -d " + i
#         subprocess.Popen(mega_align_args, shell=True).wait()
#         d = open(done,"a")
#         d.write(pmaa + '\n')
#         d.close()


infiles = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/*.fasta")

done = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/sel_done.txt"

#os.mkdir("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/alignments")
#os.mkdir("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/pos_selection")
for i in infiles:
    d = open(done,"r").readlines()
    d = [x.strip('\n') for x in d]
    print float(len(d)*100/len(infiles))
    pmaa = i.split('/')[-1][:11]
    if pmaa not in d:
        geneno = open(i,'r').read()
        if geneno.count(">") > 1:
            mega_align_args = "/usr/local/bin/megacc -a /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/clustal_align_coding.mao -o /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/alignments -s -d " + i
            subprocess.Popen(mega_align_args, shell=True).wait()
            inp = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/alignments/"+pmaa+"*.meg")[0]
            mega_ztest_args = "/usr/local/bin/megacc -a /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/zTest_coding.mao -o /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/pos_selection -d " + inp
            subprocess.Popen(mega_ztest_args, shell=True).wait()
            d = open(done,"a")
            d.write(pmaa + '\n')
            d.close()
