__author__ = 'mjohnpayne'


## counts of genes with x characteristic that are gained/lost in tm branch

## characteristics of lost gene are determined by majority type for ortho_group

## - 1 get genes into dictionaries with orthogroup and lists of characteristics

## - 2 get most common characteristic for each group

## - 3 go to badirate output and sum up gains and losses for each characteristic

## - 4 make graph of gains/losses of characteristics for each type (i.e. merops, p450, dbcan, interpro etc)


import glob
from time import sleep as sl
from collections import Counter
import seaborn as sns
from matplotlib import pyplot as plt


def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]


annot_lists = glob.glob("/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_functional_annot/*_descriptions.txt")

bdrate_changes = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/eurot_gene_gain_loss_from_badirate.txt","r").readlines()
num_change = {}

outfile = open("/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_functional_annot/orthogroup_gain_loss_functional.txt","w")

for i in bdrate_changes[1:]:
    col = i.split('\t')
    num_change[col[0]] = int(col[3])


dict = {}

for i in annot_lists:
    spec = i.split('/')[-1][:-17]
    print spec
    info = open(i,"r").read().split('\r')
    c = 1
    for line in info[1:]:
        col = line.split('\t')
        shl2 = "inconclusive"
        try:
            if col[11] == "#N/A":
                continue
            elif float(col[11]) > 0.4:
                shl2 = col[10]
        except:
            print col
        if col[13] == "-":
            dict[spec+"_singleton_" + str(c)] = {"dbCAN":[col[3]],"Merops":[col[5]],"SignalP":[col[8]],"TMHMM":[col[9]],"Sherloc2":[shl2],"p450":[col[12]],"transportTP":[col[14]],"Smurf_bb":[col[17]]}
            c +=1
        else:
            if col[13] not in dict:
                dict[col[13]] = {"dbCAN":[col[3]],"Merops":[col[5]],"SignalP":[col[8]],"TMHMM":[col[9]],"Sherloc2":[shl2],"p450":[col[12]],"transportTP":[col[14]],"Smurf_bb":[col[17]]}
            else:
                dict[col[13]]["dbCAN"].append(col[3])
                dict[col[13]]["Merops"].append(col[5])
                dict[col[13]]["SignalP"].append(col[8])
                dict[col[13]]["TMHMM"].append(col[9])
                dict[col[13]]["Sherloc2"].append(shl2)
                dict[col[13]]["p450"].append(col[12])
                dict[col[13]]["transportTP"].append(col[14])
                dict[col[13]]["Smurf_bb"].append(col[17])

annots = {}
for i in dict:
    annots[i] = {}
    for j in dict[i]:
        annots[i][j] = Most_Common(dict[i][j])

changes = {"dbCAN":{},"Merops":{},"SignalP":{},"TMHMM":{},"Sherloc2":{},"p450":{},"transportTP":{},"Smurf_bb":{}}

for i in annots:
    diff = 0
    if i in num_change:
        diff = num_change[i]
    elif "T_marneffei_singleton" in i:
        diff = 1
    for j in annots[i]:
        ##Use dbcan to get higher level category counts

        # if j == "dbCAN":
        #     type = ''.join([x for x in annots[i][j] if not x.isdigit()])

        ##Use Merops to get higher level category counts

        # if j == "Merops":
        #     type = annots[i][j][0]


## change 2nd if to "p450" for binary p450 output
        if j == "p450":
            if annots[i][j] != "-":
                type = annots[i][j]
            else:
                type = annots[i][j]
        else:
            type = annots[i][j]
        if type not in changes[j]:
            changes[j][type] = [0,0]
            if diff > 0:
                changes[j][type][0] += diff
            elif diff < 0:
                changes[j][type][1] += diff
        else:
            if diff > 0:
                changes[j][type][0] += diff
            elif diff < 0:
                changes[j][type][1] += diff

for i in changes:
    if "-" in changes[i]:
        del(changes[i]["-"])

for i in changes:
    pos = [changes[i][j][0] for j in changes[i].keys()]
    neg = [changes[i][j][1] for j in changes[i].keys()]
    cats = [j for j in changes[i].keys()]
    outfile.write(i+"\t\t\n")
    for j in changes[i]:
        cat = j
        pos = str(changes[i][j][0])
        neg = str(changes[i][j][1])
        outfile.write(cat+"\t"+pos+"\t"+neg+"\n")
    outfile.write("\n")
outfile.close()



