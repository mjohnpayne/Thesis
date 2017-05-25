__author__ = 'mjohnpayne'

import glob
import subprocess


### run delly on all bam files in list for each mut type in typels

inbams = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/CI_trimmed_fastqs/*_sort.bam")





vcf_header = """##fileformat=VCFv4.1
##fileDate=20151118
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="PE support below 3 or mapping quality below 20.">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=CONTROL,Number=1,Type=Integer,Description="Control variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled genotype likelihoods for RR,RA,AA genotypes">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Raw high-quality read counts for the SV">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\n"""


def run_delly(bams):
    typels = ["DEL","INV","DUP","TRA"]
    for i in bams:
        strain = i.split('/')[-1].replace("_sort.bam","")
        for j in typels:
            delly_args = "/usr/local/bin/delly -o /Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/"+strain +"_"+j+".vcf -g /Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta -t "+ j +" "+ i
            subprocess.Popen(delly_args, shell=True).wait()

def read_depths():
    inf = open('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/depthofcoverage_stats/CI_median_readdepths.txt',"r")
    rd = {}
    for i in inf:
        col = i.strip('\n').split('\t')
        rd[col[0]] = float(col[1])
    return rd

## filters delly predictions:
#must pass bot genotyping and basic filter
#must have split read support and have more than 20 split read supports reads
## overall size must be < 500Kb

def filt_in(indelly):
    outls = []
    strain = indelly.split('/')[-1][:-8]
    rds = read_depths()
    inf = open(indelly,"r")
    for i in inf:
        if i[0] != "#":
            col = i.split("\t")
            details = col[7].split(';')
            #print details
            pe_reads = int(details[8].replace('PE=',''))
            st = int(col[1])
            end = int(details[6][4:])
            endcont = details[5][5:]
            pelim = 40
            if strain in rds:
                pelim = rds[strain]*0.35
            if col[6] == "PASS" and details[0] == "PRECISE" and pe_reads > pelim and ":LowQual:" not in i:
                sup_reads = int(details[10].replace('SR=',''))
                if sup_reads >= 20:
                    if "TRA" not in col[4]:
                        if col[0] == endcont and end-st < 500000:
                            outls.append(i)
                    else:
                        outls.append(i)

    return outls

def remove_dupes(vcflis):
    outlis = []
    poslis = {}
    for i in vcflis:
        pres = 'N'
        col = i.split('\t')
        cont = col[0]
        pos = int(col[1])
        if cont in poslis:
            for j in poslis[cont]:
                j = float(j)
                if j-30 < pos < j+30:
                    pres = 'Y'
        if pres == 'N':
            if cont not in poslis:
                poslis[cont] = [pos]
            else:
                poslis[cont].append(pos)
            outlis.append(i)
    return outlis

def make_rep_dict(reps):
    repeats = {}
    for i in reps:
        if i[0] != "#":
            col = i.split('\t')
            cont = col[0][:-2]
            if cont not in repeats:
                repeats[cont] = [[int(col[3]),int(col[4])]]
            else:
                repeats[cont].append([int(col[3]),int(col[4])])
    return repeats

## removes delly annotations that have break points in annotated repetitive DNA

def remove_overlap_reps(rep,vcfls):
    outls = []
    for i in vcfls:
        col = i.split('\t')
        cont = col[0]
        pos = int(col[1])
        details = col[7].split(';')
        end = int(details[6][4:])
        endcont = details[5][5:]
        repetitive = 'N'
        for j in rep[cont]:
            if j[0] < pos < j[1]:
                repetitive = 'Y'
        for j in rep[endcont]:
            if j[0] < end < j[1]:
                repetitive = 'Y'
        if repetitive == 'N':
            outls.append(i)
    return outls

def merge_types(vars):
    repeats = open("/Volumes/MP_HD/repDNA_data/repeatmasker/Pm_all_genome.fasta.out.gff","r")
    reps = make_rep_dict(repeats)
    repeats.close()
    typels = ["DEL","INV","DUP","TRA"]
    for j in typels:
        filtls = []
        for i in vars:
            if j in i:
                filtls += filt_in(i)
        filtls = remove_dupes(filtls)
        filtls = remove_overlap_reps(reps,filtls)
        outfile = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filt_precise_all_" + j + ".vcf","w")
        outfile.write(vcf_header + "".join(filtls))

def filter_individually(vars):
    outf = vars.split("/")
    o1 = outf[-1][:-4]
    o2 = '/'.join(outf[:-2]) + "/filtered_vcf/filtered_" + o1 + ".vcf"
    repeats = open("/Volumes/MP_HD/repDNA_data/repeatmasker/Pm_all_genome.fasta.out.gff","r")
    reps = make_rep_dict(repeats)
    repeats.close()
    filtls = filt_in(vars)
    filtls = remove_dupes(filtls)
    filtls = remove_overlap_reps(reps,filtls)
    outfile = open(o2,"w")
    outfile.write(vcf_header + "".join(filtls))



#run_delly(inbams)

allvar = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/raw_vcf/*.vcf")

#merge_types(allvar)
for i in allvar:
    filter_individually(i)