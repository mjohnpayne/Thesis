__author__ = 'mjohnpayne'

import glob

vcfls = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/*.vcf")

def parse_gff():
    gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3",'r')
    genes = {}
    for line in gff:
        if line[0] != "#":
            col = line.split("\t")
            if col[2] == 'gene':
                st = int(col[3])
                en = int(col[4])
                orient = col[6]
                id = col[8].split(";")[1].replace("Name=","")
                if col[0] not in genes:
                    genes[str(col[0])] = [(st,en,orient,id)]
                else:
                    genes[str(col[0])].append((st,en,orient,id))
    return genes


def process_dels(dst,den,gst,gen,gorient):
    delperc = 0.0
    delpart = ""
    if dst < gen < den and gst < dst:
        if gorient == "+":
            lendel = gen-dst
            genlen = gen-gst
            delperc = (float(lendel)/genlen)*100
            delpart = "3p"
        if gorient == "-":
            lendel = gen-dst
            genlen = gen-gst
            delperc = (float(lendel)/genlen)*100
            delpart = "5p"
    elif dst < gst < den and gen > den:
        lendel = den-gst
        genlen = gen-gst
        delperc = (float(lendel)/genlen)*100
        if gorient == "+":
            delpart = "5p"
        elif gorient == "-":
            delpart = "3p"
    elif dst < gst  and den > gen:
        delperc = 100.0
        delpart = "all"
    elif dst > gst and den < gen:
        lendel = den-dst
        genlen = gen-gst
        delperc = (float(lendel)/genlen)*100
        delpart = "internal"
    return delperc,delpart

# def process_inv(ist,ien,gst,gen,gorient):
#     invperc = 0.0
#     invpart = ""
#     if ist < gen < ien and gst < ist:
#         if gorient == "+":
#             leninv = gen-ist
#             genlen = gen-gst
#             invperc = (float(leninv)/genlen)*100
#             invpart = "3p"
#         if gorient == "-":
#             leninv = gen-ist
#             genlen = gen-gst
#             invperc = (float(leninv)/genlen)*100
#             invpart = "5p"
#     elif ist < gst < ien and gen > ien:
#         leninv = ien-gst
#         genlen = gen-gst
#         invperc = (float(leninv)/genlen)*100
#         if gorient == "+":
#             invpart = "5p"
#         elif gorient == "-":
#             invpart = "3p"
#     elif ist > gst and ien < gen:
#         leninv = ien-ist
#         genlen = gen-gst
#         invperc = (float(leninv)/genlen)*100
#         invpart = "internal"
#     return invperc,invpart

def process_tra(tpoint,gst,gen,gorient):
    lossperc = 0.0
    if gorient == "+":
        losslen = gen-tpoint
        glen = gen-gst
        lossperc = (float(losslen)/glen)*100
    elif gorient == "-":
        losslen = tpoint-gst
        glen = gen-gst
        lossperc = (float(losslen)/glen)*100
    return lossperc,"N/A"



def parse_vcf(invcf,genepos):
    inf = open(invcf,"r")
    muts = {}
    for i in inf:
        if i[0] != "#":
            col = i.split("\t")
            type = col[4][1:-1]
            details = col[7].split(';')
            cont = str(col[0])
            pe_reads = int(details[8].replace('PE=',''))
            st = int(col[1])
            end = int(details[6][4:])
            endcont = details[5][5:]
            ## if start > inv start or if start <del end
            if type == "DEL":
                for j in genepos[cont]:
                    if st < j[0] <end or st < j[1] < end:
                        perc_del,del_part = process_dels(st,end,j[0],j[1],j[2])
                        if j[3] not in muts:
                            muts[j[3]] = [("Deletion",perc_del,del_part)]
                        else:
                            muts[j[3]].append(("Deletion",perc_del,del_part))
            if type == "INV":
                for j in genepos[cont]:
                    if j[0] < st < j[1]:
                        perc_inv,inv_part = process_tra(st,j[0],j[1],j[2])
                        if j[3] not in muts:
                            muts[j[3]] = [("inversion",perc_inv,inv_part)]
                        else:
                            muts[j[3]].append(("inversion",perc_inv,inv_part))
                for j in genepos[endcont]:
                    if j[0] < end < j[1]:
                        lperc,ltype = process_tra(end,j[0],j[1],j[2])
                        if j[3] not in muts:
                            muts[j[3]] = [("inversion",lperc,ltype)]
                        else:
                            muts[j[3]].append(("inversion",lperc,ltype))
            if type == "TRA":
                for j in genepos[cont]:
                    if j[0] < st < j[1]:
                        lperc,ltype = process_tra(st,j[0],j[1],j[2])
                        if j[3] not in muts:
                            muts[j[3]] = [("Translocation",lperc,ltype)]
                        else:
                            muts[j[3]].append(("Translocation",lperc,ltype))
                for j in genepos[endcont]:
                    if j[0] < end < j[1]:
                        lperc,ltype = process_tra(end,j[0],j[1],j[2])
                        if j[3] not in muts:
                            muts[j[3]] = [("Translocation",lperc,ltype)]
                        else:
                            muts[j[3]].append(("Translocation",lperc,ltype))
            if type == "DUP":
                for j in genepos[cont]:
                    if st < j[0] and end > j[1]:
                        if j[3] not in muts:
                            muts[j[3]] = [("Duplication",100.0,'N/A')]
                        else:
                            muts[j[3]].append(("Duplication",100.0,'N/A'))

    return muts

def outsum_eachstrain(vcfs,type,of):
    muts = {}
    outf = open(of,"w")
    outf.write("GeneID")
    slis = []
    for i in vcfs:
        if type in i:
            strain = i.split('/')[-1][9:-8]
            slis.append(strain)
            #outf.write("\t"+strain)
            if strain not in muts:
                muts[strain] = {}
                for t in gdict:
                    for j in gdict[t]:
                        muts[strain][j[3]] = []
            mutats = parse_vcf(i,gdict)
            for j in mutats:
                muts[strain][j]=mutats[j]
    outf.write('\t'+'\t'.join(sorted(slis))+'\n')
    for t in gdict:
        for j in gdict[t]:
            #print j[3]
            outf.write(j[3])
            for x in sorted(slis):
                if len(muts[x][j[3]]) > 0:
                    #print muts[x][j[3]]
                    outf.write("\t" + str(muts[x][j[3]][0][1])+","+muts[x][j[3]][0][2])
                else:
                    outf.write("\t-")
            outf.write("\n")

gdict = parse_gff()

#mutats = parse_vcf("/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_precise_alltypes.vcf",gdict)

outsum_eachstrain(vcfls,"INV","/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/filtered_INV_summary.txt")












