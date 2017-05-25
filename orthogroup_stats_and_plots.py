__author__ = 'mjohnpayne'

from numpy import average as av
from numpy import std
from numpy import sqrt
from time import sleep as sl
import seaborn as sns
from matplotlib import pyplot as plt

orthogroups = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/orthomcl_ortho_groups.txt'

#of = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/cluster_class_data.txt','w')

gffs = {'Pm':'/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_databases/pm_wt_dbs/pmfa1_working_models_fix.gff','Ts':'/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_databases/ts_wt_dbs/tsta1_working_models.gff3',"Tf":"/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_databases/tf_wt_dbs/TF_vel_pfams_para_genome_fix_rename.gff","Pf":"/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_databases/pf_wt_dbs/Pf_vel_denovo_rename.gff"}

group_cats = ['Absolutely conserved','Conserved','PM_only','TF_only','PF_only','TS_only','Other']


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


def get_orthoids(og,type):
    orthos = []
    allgenes = []
    og = open(og,'r')
    for line in og:
        id = ''
        genes = []

        col = line.strip('\n').split(' ')
        id = col[0].replace(':','')
        genes = col[1:]
        ngenes = []
        for i in genes:
            i = rn(i)
            ngenes.append(i)
        if count_ortho_nos(line) == type:
            orthos += ngenes
        allgenes += ngenes
    og.close()
    return orthos, allgenes


def count_ortho_nos(ln):
    if ln.count('Pma|') == 1 and ln.count('Tst|') == 1 and ln.count('Tfl|') == 1 and ln.count('Pfu|') == 1:
        return 'Absolutely conserved'
    elif ln.count('Pma|') >= 1 and ln.count('Tst|') >= 1 and ln.count('Tfl|') >= 1 and ln.count('Pfu|') >= 1:
        return 'Conserved'
    elif ln.count('Pma|') >= 1 and ln.count('Tst|') == 0 and ln.count('Tfl|') == 0 and ln.count('Pfu|') == 0:
        return 'PM_only'
    elif ln.count('Pma|') == 0 and ln.count('Tst|') >= 1 and ln.count('Tfl|') == 0 and ln.count('Pfu|') == 0:
        return 'TS_only'
    elif ln.count('Pma|') == 0 and ln.count('Tst|') == 0 and ln.count('Tfl|') >= 1 and ln.count('Pfu|') == 0:
        return 'TF_only'
    elif ln.count('Pma|') == 0 and ln.count('Tst|') == 0 and ln.count('Tfl|') == 0 and ln.count('Pfu|') >= 1:
        return 'PF_only'
    else:
        return 'Other'


def get_gene_conigs(gfin):
    gf = open(gfin,'r')
    conts = {}
    gen = {}
    for line in gf:
        if not line.startswith('#'):
            col = line.strip('\n').split('\t')
            if col[2] == 'gene':
                det = col[8].split(';')
                pmaa=''
                if col[8].startswith('ID=gene'):
                    pmaa = det[1][5:]
                else:
                    pmaa = det[0][3:]
                st = int(col[3])
                en = int(col[4])
                orient = col[6]
                cont = ''
                if 'TF_vel' in gfin:
                    cont = 'Tf' + col[0]
                else:
                    cont = col[0]
                if cont in conts:
                        conts[cont] += [pmaa]
                elif cont not in conts:
                        conts[cont] = [pmaa]
                gen[pmaa] = [st,en,orient]
    for cont in conts:
        conts[cont] = sorted(conts[cont])
    gf.close()
    gf = open(gfin,'r').read().split('\tgene\t')
    for i in gf:
        if 'CDS' in i:
            cdss = i.count('\tCDS\t')
            pmaa = ''
            i = i.replace('\n','\t')
            col = i.split('\t')
            det = col[5].split(';')
            if col[5].startswith('ID=gene'):
                pmaa = det[1][5:]
            else:
                pmaa = det[0][3:]
            gen[pmaa].append(cdss)

    return conts,gen


def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

def get_stats(gendict,genlis,spec,type,algen):
    genlen = []
    cdsno = []
    for i in genlis:
        if spec in i:
            if i in gendict:
                gene = gendict[i]
                if gene[1]-gene[0] < 0:
                    genlen.append(0)
                else:
                    genlen.append(gene[1]-gene[0])
                cdsno.append(gene[3])
    if 'only' in type:
        for i in gendict:
            if spec in i and type[:2] in i:
                if i not in algen:
                    # if spec == 'PM':
                    #     print i
                    gene = gendict[i]
                    if gene[1]-gene[0] < 0:
                        genlen.append(0)
                    else:
                        genlen.append(gene[1]-gene[0])
                    cdsno.append(gene[3])
    if len(genlen) > 0:
        return genlen, cdsno
    else:
        return 'N/A', 'N/A'


def stats(og,gfls,type,spec):
    contigs = {}
    genes = {}
    lendict = {}
    exondict = {}
    for gff in gfls:
        c2,g2 = get_gene_conigs(gffs[gff])
        contigs = merge_two_dicts(contigs,c2)
        genes = merge_two_dicts(genes,g2)
    geneids,allgen = get_orthoids(og,type)
    lengths,cdsnum = get_stats(genes,geneids,spec,type,allgen)
    return lengths,cdsnum

species = ['PM','TF','PF','TS']

#of.write('Group type\tspecies\tgene number\taverage length\tstdev\tsem\texon number average\texon number stdev\texon number sem\n')


lendict = {}
exondict = {}

lentypes = {}
cdstypes = {}
speclist = {}
lentypes["Species specific"] = []
cdstypes["Species specific"] = []
speclist["Species specific"] = []
for i in group_cats:
    print i
    lendict[i] = {}
    exondict[i] = {}
    lenlist = []
    cdslist = []
    if "only" not in i:
        lentypes[i] = []
        cdstypes[i] = []
        speclist[i] = []
    for j in species:
        length,cds = stats(orthogroups,gffs,i,j)
        lendict[i][j] = length
        exondict[i][j] = cds
        print len(length)
        if "only" in i and len(length) > 10:
            speclist["Species specific"] += len(length)*[j]
            lentypes["Species specific"] += length
            cdstypes["Species specific"] += cds
        elif "only" in i and len(length) < 10:
            continue
        else:
            speclist[i] += len(length)*[j]
            lentypes[i] += length
            cdstypes[i] += cds
traces = []



fig, ax = plt.subplots()
order = ["Absolutely conserved","Conserved","Species specific",'Other']
values = []
species = []
cat = []
for i in order:
    values += lentypes[i]
    cat += len(lentypes[i])*[i]
    species += speclist[i]

ax = sns.boxplot(x=cat, y=values, hue=species, palette="Set1")
plt.ylim([0,12500])

fig.show()
#fig.savefig('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/cluster_class_gene_lengths.pdf')


# fig, ax = plt.subplots()
# order = ["Absolutely conserved","Conserved","Species specific",'Other']
# values = []
# species = []
# cat = []
# for i in order:
#     values += cdstypes[i]
#     cat += len(cdstypes[i])*[i]
#     species += speclist[i]
#
# ax = sns.boxplot(x=cat, y=values, hue=species, palette="Set1")
# plt.ylim([0,20])
#
# fig.savefig('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/cluster_class_gene_exons.pdf')