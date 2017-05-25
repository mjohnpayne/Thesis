__author__ = 'mjohnpayne'

from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace
import math
from time import sleep as sl

inbd = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/eurot_gain_loss_nofam_nosingle.bd",'r')
outfile = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/eurot_gene_gain_loss_from_badirate.txt",'w')

outfile.write("Group\tNo at immediate branch\tNo in Pm\tChange from branch to Pm\n")

for line in inbd:
    if "		eurot" in line:
        pmno = line[line.find("Pma_")+4:]
        bno = pmno[pmno.find(")")+1:]
        pmno = float(pmno[:pmno.find(":")])
        bno = float(bno[:bno.find(":")])
        change_frm_branch = bno - pmno
        group = line[line.find("eurot_group"):line.find("(")-1]
        outfile.write(group+'\t'+str(int(bno)) + '\t' + str(int(pmno)) + '\t' + str(int(pmno-bno)) + '\n')
outfile.close()