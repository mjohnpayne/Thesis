__author__ = 'mjohnpayne'

import glob
from time import sleep as sl



inlist = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/TCDB/transportTP/*PredictedTransporters.tab")

outf = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/TCDB/transportTP/"

def parse_tmhmmout(inls,outpref):
    for i in inls:
        inf = open(i,"r").readlines()
        spec = i.split('/')[-1][:2]
        out = open(outf+spec+"transportTP_parsed.txt","w")
        out.write('GeneID\ttransporter id\ntransporter description\n')
        for j in inf[16:]:
            col = j.split('\t')
            pmaa = col[0]
            famname = col[4]
            famdesc = col[9]
            out.write(pmaa+'\t'+famname+'\t'+famdesc+'\n')
        out.close()


parse_tmhmmout(inlist,outf)