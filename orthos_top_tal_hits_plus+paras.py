__author__ = 'mjohnpayne'
import subprocess
from Bio import SeqIO
import os
import glob
import re
from time import sleep as sl

orthogroup_fastas = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/4species_orthogroups_fastas/all_groups/*.fasta")

orthos_by_gene = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_groups_by_gene.txt','r')

all_talaro = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/all_talaro_proteins.fasta"

all_talaro_cds = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/all_talaro_cds.fasta"

done = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/fasta_done.txt"

all = SeqIO.parse(all_talaro,"fasta")
all_seqs = {}

for i in all:
    all_seqs[i.id] = i

allc = SeqIO.parse(all_talaro_cds,"fasta")
all_cds_seqs = {}
counts = {}
for i in allc:
    all_cds_seqs[i.id] = i
#     if i.id not in counts:
#         counts[i.id] = 1
#     else:
#         counts[i.id] += 1
# c = 0
# for i in counts:
#     if counts[i] > 1 and "TFLA" in i:
#         print i
#         c +=1
# print c


def rn(gene):
    if 'TF' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    elif 'Pf' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'PFUN_' + no*'0' + gene[3:] + '0'
    else:
        gene = gene[4:]
    return gene






def make_blastdb(infasta):
    mkdbargs = "/usr/bin/makeblastdb -in " + infasta + " -dbtype prot"
    subprocess.Popen(mkdbargs, shell=True).wait()

def run_blast_return_top(db,fasta):
    blast_args = '/opt/local/bin/blast -d ' + db + ' -i ' + fasta + ' -m 8 -p blastp'

    blast_out = subprocess.Popen(blast_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

    out = blast_out.communicate()[0]
    out = str(out)
    out = out.split('\n')
    outls = []
    for line in out[1:-1]:
        line = line.split('\t')
        outls.append(rn(line[1]))
    idstring = "".join(outls)
    if "TFLA" in idstring:
        for i in outls:
            if "TFLA" in i:
                return i
    elif "PFUN" in idstring:
        for i in outls:
            if "PFUN" in i:
                return i
    elif "TSTA" in idstring:
        for i in outls:
            if "TSTA" in i:
                return i
    else:
        return "NONE"


def run_blast_return_top_talaros(db,fasta):
    blast_args = '/opt/local/bin/blast -d ' + db + ' -i ' + fasta + ' -m 8 -p blastp'

    blast_out = subprocess.Popen(blast_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

    out = blast_out.communicate()[0]
    out = str(out)
    out = out.split('\n')
    outls = []
    for line in out[1:-1]:
        line = line.split('\t')
        outls.append(line[1])
    idstring = "".join(outls)
    tf = idstring.find("TFLA")
    tf = idstring[tf:tf+11]
    pf = idstring.find("PFUN")
    pf = idstring[pf:pf+11]
    ts = idstring.find("TSTA")
    ts = idstring[ts:ts+11]
    return [tf,pf,ts]



##return fasta with 1st gene of interest, followed by closest talaro in ortho group, followed by paralogues, followed by top TF,PF,TS hits (with closest talaro removed)

def process_orthogroups(ingroup):
    group1 = SeqIO.parse(ingroup,"fasta")
    group = SeqIO.parse(ingroup,"fasta")
    ids = []
    for i in group1:
        ids += [rn(i.id)]
    make_blastdb(ingroup)
    for i in group:
        if "PMAA" in i.id:
            #paras = [str(x.id[4:]) for x in group if "PMAA" in x.id]
            paras = [x for x in ids if 'PMAA' in x]
            paras.remove(i.id[4:])
            i.id = i.id[4:]
            i.description = ""
            SeqIO.write(i,"tmp.fasta","fasta")
            top = run_blast_return_top(ingroup,"tmp.fasta")
            talaros = run_blast_return_top_talaros(all_talaro,"tmp.fasta")
            outs = []
            if top != "NONE":
                outs = [i.id,top] + paras
            else:
                outs = [i.id] + paras
            for j in talaros:
                if j != "" and j not in outs:
                    outs += [j]
            os.remove("tmp.fasta")
            SeqIO.write([all_cds_seqs[x] for x in outs],"/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_paras_and_talaros/" + i.id + "_paras_talaros.fasta","fasta")
    os.remove(ingroup + ".phr")
    os.remove(ingroup + ".pin")
    os.remove(ingroup + ".psq")
c = 0

for i in orthogroup_fastas:
    d = open(done,"r").readlines()
    d = [x.strip('\n') for x in d]
    if i not in d:
        process_orthogroups(i)
        print c*100/len(orthogroup_fastas),"%"
        d = open(done,"a")
        d.write(i + '\n')
        d.close()
    c +=1
#top = run_blast_return_top("/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/pm_proteins_with_byss.fasta","/Users/mjohnpayne/Documents/PhD/ASPS/Asp_blast_fro_tree/PMAA_090410.fasta")

