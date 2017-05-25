__author__ = 'mjohnpayne'

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from time import sleep as sl
import re


ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')

#ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/tsta1_working_models.gff3','r')

#record = SeqIO.read("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pm_genbank_gene_containing_only.gb", "genbank")

#helicases = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/recQ_helicases/recQ_helicase_ts_gene_ids.txt",'r').read().split('\r')

helicases = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/recQ_helicases/recQ_helicase_gene_ids.txt",'r').read().split('\r')
# pops = open("/Users/mjohnpayne/Documents/PhD/ASPS/pop_ids.txt",'r').read().split('\n')
# byss = open('/Users/mjohnpayne/Documents/PhD/bys/complete_list_bys_ids_f.txt','r').read().split('\n')

#print helicases
#print pops

genome = {}
for i in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
    genome[i.id] = i.seq
size = 40000
sname = str(float(size)/1000)
ends5 = []
ends3 = []
for i in genome:
    l = len(genome[i])
    if l > 50000:
        ends5 += [(i,0,size)]
        ends3 += [(i,l-size,l)]

print ends3
contigs = {}
gene = {}
for line in ingff:
        col = line.strip('\n').split('\t')
        if 'ID=gene' in line:
            det = col[8].split(';')
            pmaa = det[1][5:]
            st = int(col[3])
            en = int(col[4])
            orient = col[6]
            cont = col[0]
            if pmaa not in gene:
                if cont in contigs:
                        contigs[cont] += [pmaa]
                elif cont not in contigs:
                        contigs[cont] = [pmaa]

                gene[pmaa] = [cont,st,en,orient]

tel_rep = ['TTAGGG','CCCTAA','TTAGGA','TCCTAA'] ## -ve

TTAGGG = []
TTAGGA = []

for i in tel_rep:
    for pos in ends5:
        telomere_sites = []
        telomere_sites = [m.start() for m in re.finditer(i,str(genome[pos[0]][pos[1]:pos[2]]))]
        for j in telomere_sites:
            j += pos[1]
            name = i+"_"+pos[0]+"_"+str(j)
            if i == 'TTAGGG':
                TTAGGG.append(name)
                gene[name] = [pos[0],j,j+6,"+"]
            elif i == 'CCCTAA':
                TTAGGG.append(name)
                gene[name] = [pos[0],j,j+6,"-"]
            elif i == "TTAGGA":
                TTAGGA.append(name)
                gene[name] = [pos[0],j,j+6,"+"]
            elif i == 'TCCTAA':
                TTAGGA.append(name)
                gene[name] = [pos[0],j,j+6,"-"]

for i in tel_rep:
    for pos in ends3:
        telomere_sites = []
        telomere_sites = [m.start() for m in re.finditer(i,str(genome[pos[0]][pos[1]:pos[2]]))]
        for j in telomere_sites:
            j += pos[1]
            name = i+"_"+pos[0]+"_"+str(j)
            if i == 'TTAGGG':
                TTAGGG.append(name)
                gene[name] = [pos[0],j,j+6,"+"]
            elif i == 'CCCTAA':
                TTAGGG.append(name)
                gene[name] = [pos[0],j,j+6,"-"]
            elif i == "TTAGGA":
                TTAGGA.append(name)
                gene[name] = [pos[0],j,j+6,"+"]
            elif i == 'TCCTAA':
                TTAGGA.append(name)
                gene[name] = [pos[0],j,j+6,"-"]


#pos = ('93',0,40000)




lab = False
lsize = 5
athick = 1


gdd = GenomeDiagram.Diagram("5_ends")
c = 1
for pos in ends5:
    gd_track_for_features = gdd.new_track(c, name=pos[0], greytrack=True)
    gd_feature_set = gd_track_for_features.new_set()
    for i in gene:
        if gene[i][0] == pos[0] and gene[i][1] < pos[2] and i in helicases:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=i, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=i, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        elif gene[i][0] == pos[0] and gene[i][1] < pos[2] and i in TTAGGG:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        elif gene[i][0] == pos[0] and gene[i][1] < pos[2] and i in TTAGGA:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="blue", arrowshaft_height=athick, label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="blue", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        elif gene[i][0] == pos[0] and gene[i][1] < pos[2]:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=150)

    # if pos[1] == 0:
    #     feature = SeqFeature(FeatureLocation(0, 200), strand=0)
    #     gd_feature_set.add_feature(feature, name="Scaffold end", label=lab, sigil="BOX", color="red", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
    # else:
    #     feature = SeqFeature(FeatureLocation(pos[2]-200, pos[2]), strand=0)
    #     gd_feature_set.add_feature(feature, name="Scaffold end", label=lab, sigil="BOX", color="red", arrowshaft_height=athick,  label_size = lsize, label_angle=30)

    # tel_rep = 'TAAGGG'## or in nidulans 'TAAGGG'
    # telomere_sites = [m.start() for m in re.finditer(tel_rep,str(genome[pos[0]][pos[1]:pos[2]]))]
    # for i in telomere_sites:
    #     feature = SeqFeature(FeatureLocation(i, i+5))
    #     gd_feature_set.add_feature(feature, color="blue", name='telomeric_rep',label=True, label_size = lsize)

    c +=1


gdd.draw(format='linear', pagesize='A4', fragments=1, start=0, end=size,track_size=0.5)

gdd.write("/Volumes/MP_HD/repDNA_data/telomere_figures_frm_script/5_prime_"+sname+"kb_tel_repeats_labels.pdf", "pdf")


lab = False
lsize = 5
athick = 1


gdd = GenomeDiagram.Diagram("3_ends")
c = 1
for pos in ends3:
    print pos
    gd_track_for_features = gdd.new_track(c, name=pos[0], greytrack=True)
    gd_feature_set = gd_track_for_features.new_set()
    for i in gene:
        gstart = gene[i][1] - pos[1]
        gend = gene[i][2] - pos[1]
        if gene[i][0] == pos[0] and gend > 0 and i in helicases:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=i, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=i, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        elif gene[i][0] == pos[0] and gend > 0  and i in TTAGGG:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        elif gene[i][0] == pos[0] and gend > 0  and i in TTAGGA:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="blue", arrowshaft_height=athick, label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="blue", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        elif gene[i][0] == pos[0] and gend > 0 :
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gstart, gend), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=150)

    # if pos[1] == 0:
    #     feature = SeqFeature(FeatureLocation(0, 200), strand=0)
    #     gd_feature_set.add_feature(feature, name="Scaffold end", label=lab, sigil="BOX", color="red", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
    # else:
    #     feature = SeqFeature(FeatureLocation(pos[2]-200, pos[2]), strand=0)
    #     gd_feature_set.add_feature(feature, name="Scaffold end", label=lab, sigil="BOX", color="red", arrowshaft_height=athick,  label_size = lsize, label_angle=30)

    # tel_rep = 'TAAGGG'## or in nidulans 'TAAGGG'
    # telomere_sites = [m.start() for m in re.finditer(tel_rep,str(genome[pos[0]][pos[1]:pos[2]]))]
    # for i in telomere_sites:
    #     feature = SeqFeature(FeatureLocation(i, i+5))
    #     gd_feature_set.add_feature(feature, color="blue", name='telomeric_rep',label=True, label_size = lsize)

    c +=1


gdd.draw(format='linear', pagesize='A4', fragments=1, start=0, end=size,track_size=0.5)

gdd.write("/Volumes/MP_HD/repDNA_data/telomere_figures_frm_script/3_prime_"+sname+"kb_tel_repeats_labels.pdf", "pdf")





