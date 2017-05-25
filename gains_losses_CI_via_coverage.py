__author__ = 'mjohnpayne'

from time import sleep as sl

import numpy as np

def iround(x):
    """iround(number) -> integer
    Round a number to the nearest integer."""
    return int(round(x) - .5) + (x > 0)

def gene_int(array,cop):
    return np.round(array/cop)



incov = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene coverage/cov_comparison_between_strains_corrected_for_tot_reads_for_python.txt",'r').read().split('\r')
outfile = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/gene_counts_CIs.tsv",'w')




strains = incov[0].split('\t')[1:]
#print strains

outfile.write("ID\tFRR2161\t" + "\t".join(strains) + '\n')

for i in incov[1:]:
    gene  = i.split('\t')[0]
    cov = map(float,i.split('\t')[1:])
    cov = np.asarray(cov)
    copy1 = cov[(cov>60)&(cov<160)]
    copies = []
    if len(copy1) == 0:
        arr_avg = np.average(cov)
        if arr_avg > 60:
            copies = gene_int(cov,100)
        elif arr_avg < 60:
            copies = [0]*20
    else:
        c1 = np.average(copy1)
        copies = gene_int(cov,c1)
    copies = map(str,map(int,list(copies)))
    outfile.write(gene + '\t1\t'+ '\t'.join(copies) + "\n")

outfile.close()


# ((((4059:2.3550391395962844E-5,3840:1.9482814206100593E-5)10:6.128652506132148E-5,F4:1.3235217487967088E-4)10:2.064778428707447E-4,3841:2.9743397381548747E-4)10:5.487300437392128E-5,(((((BR2SD2:2.789478316937542E-5,BR2:2.6331002328641466E-5)4:4.482127137874544E-7,BR2SD1:2.480514311379331E-5)10:1.8944295611177786E-4,'043':2.6022292104017686E-4)10:8.903564191307955E-5,'702':2.9958960308882276E-4)10:1.3421716924330603E-5,((((('203SD4':2.5734406151208878E-5,'203SD3':2.6578711455707392E-5)9:1.2336708353447457E-6,'203':3.198915760231304E-5)10:1.9206559436700217E-6,'203SD1':3.268358909113229E-5)10:1.9194419426567596E-4,(('027SD2':2.6113249297561895E-5,'027SD1':2.6660324948421902E-5)10:8.113617113169453E-6,'027':3.368092194938595E-5)10:1.960869077183465E-4)10:5.5059292991640324E-5,((('3871':1.702603925147246E-4,'3482':1.4579898867410893E-4)10:6.075399178693359E-5,HR2:2.1002943347092347E-4)10:4.618659043765621E-5,'012':2.7027813566345117E-4)10:8.218538065854104E-5)10:1.844462672450518E-5)10:5.487300437392128E-5):0.0;
