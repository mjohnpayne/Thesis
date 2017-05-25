__author__ = 'mjohnpayne'

import glob
from time import sleep as sl
import re


outfile = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/pos_sel_tests.txt",'w')
meglist = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/pos_selection/*.meg")#['/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_and_top_positive_selection_tests/PMAA_002350_top_ortho-21171-40505.meg']

outfile.write("Tm gene\tClosest ortho\tdn-ds\tpvalue\n")



def megaparse(inmeg):
    #print inmeg
    tmp = open(inmeg,"r").read()
    ## pull out rows with gene and number assignment pairs. Make list of gene ids
    start = tmp.find("1] #")
    end = tmp.find("[   ")-2
    list = tmp[start:end].split('\n')
    ids = [x[-11:] for x in list]
    #identify PMAA in pos 1 and then first non pm gene following
    pmaa = ids[0]
    comp = ""
    pn = 1
    cn = 0
    for i in range(1,len(ids)):
            if 'PMAA' not in ids[i]:
                comp = ids[i]
                cn = i+1
                break
    #extract table and store in list of lists i.e. array
    table = tmp[end:].split('\n')
    table = [x[6:] for x in table]
    sep = table[2].find("2")-table[2].find("1")
    table = [[""]+[x[i:i+sep].replace("[","").replace("]","").replace(" ","") for i in range(0, len(x), sep)] for x in table]
    table = table[2:-2]
    pval = table[cn][pn]
    dnds = table[pn][cn]
    #print pmaa,comp,dnds,pval
    outfile.write(pmaa + '\t' + comp + '\t' + str(dnds) + '\t' + str(pval) + '\n')


for i in meglist:
    megaparse(i)

outfile.close()

#     print i
#     tmp = open(i,"r").readlines()
#     pval= ''
#     dnds = ''
#     comp = 1
#     for j in range(len(tmp)):
#         line = tmp[j]
#         if "[ 1] #" in line or "[1] #" in tmp[j]:
#             pmgene = line[line.find("1] #")+4:line.find("1] #")+15]
#             nline = tmp[j+1]
#             gene = nline[nline.find("] #")+3:nline.find("] #")+14]
#
#         elif "[ 1]    " in line or "[1]    " in line:
#             dnds = line[4:][line[4:].find('[')+1:line[4:].find(']')-1]
#         elif "[ 2]    " in line or "[2]    " in line:
#             pval = line[7:17]
#     print pmgene,gene,dnds,pval
#     # p1 = tmp.find("1] #")
#     # pmgene = tmp[p1+4:p1+15]
#     # p2 = tmp.find("2] #")
#     # orthgene = tmp[p2+4:p2+15]
#     # p3 = tmp.find("1]     ")
#     # dline = tmp[p3+4:p3+33]
#     # dnds = dline[dline.find('[')+1:dline.find(']')-1]
#     # p4 = tmp.find("2]   ")
#     # pline = tmp[p4:p4+33]
#     # pval = pline[7:17]
#     # #,orthgene,dnds,pval
#     sl(0.2)
# #     outfile.write(pmgene + '\t' + orthgene + '\t' + dnds + '\t' + pval + '\n')
# # outfile.close()