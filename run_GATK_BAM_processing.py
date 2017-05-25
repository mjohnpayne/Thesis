__author__ = 'mjohnpayne'

import sys
import subprocess


inbam = sys.argv[1] #'/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam/203_sort.bam'

print inbam.split('/')[-1]

print '\n\nStarting re-sort\n'

sortout = '/'.join(inbam.split('/')[:-1]) + '/GATK_processing/' + inbam.split('/')[-1][:-9] + '_re_sort.bam'


sort_args = 'java -Xmx10g -jar /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/picard-tools-1.63/picard-tools-1.63/SortSam.jar INPUT=' + inbam + ' OUTPUT=' + sortout + ' SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT'

subprocess.Popen(sort_args, shell=True).wait()



print '\n\n\n\nStarting Dedup\n'

dedup_out = sortout[:-12] + '_re_dedup_sort.bam'
dedup_met = sortout[:-12] + '_dedup_metrics.txt'

dedup_args = 'java -Xmx10g -jar /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/picard-tools-1.63/picard-tools-1.63/MarkDuplicates.jar INPUT=' + sortout + ' OUTPUT=' + dedup_out + ' METRICS_FILE=' + dedup_met + ' VALIDATION_STRINGENCY=SILENT'

subprocess.Popen(dedup_args, shell=True).wait()



print '\n\n\n\nStarting read group addition\n'

rg_out = sortout[:-12] + '_rg_dedup_re_sort.bam'
rg_args = 'java -Xmx10g -jar /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/picard-tools-1.63/picard-tools-1.63/AddOrReplaceReadGroups.jar INPUT=' + dedup_out + ' OUTPUT=' + rg_out + ' RGID=group1 RGLB= lib1 RGPL=illumina RGPU=unit1 RGSM=sample1 VALIDATION_STRINGENCY=SILENT'

subprocess.Popen(rg_args, shell=True).wait()



print '\n\n\n\nStarting bam index\n'

idx_args = 'java -Xmx10g -jar /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/picard-tools-1.63/picard-tools-1.63/BuildBamIndex.jar INPUT=' + rg_out + ' VALIDATION_STRINGENCY=SILENT'

subprocess.Popen(idx_args, shell=True).wait()



subprocess.Popen('rm ' + dedup_out, shell=True).wait()




print '\n\n\n\nStarting RealignerTargetCreator\n'

subprocess.Popen('rm ' + sortout, shell=True).wait()
intout = rg_out.replace('.bam','.intervals')
args = 'java -Xmx10g -jar /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/GATK+Queue/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T RealignerTargetCreator -fixMisencodedQuals -R /Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta -I '+ rg_out + ' -o ' + intout
subprocess.Popen(args, shell=True).wait()




print '\n\n\n\nStarting IndelRealigner\n'

realign_out = sortout[:-12] + '_re_dedup_sort_realign.bam'
args = 'java -Xmx10g -jar /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/GATK+Queue/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T IndelRealigner -fixMisencodedQuals -R /Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta -targetIntervals ' + intout + ' -I '+ rg_out + ' -o ' + realign_out

subprocess.Popen(args, shell=True).wait()




print '\n\n\n\nStarting BaseRecalibrator\n'

recal_tbl = sortout[:-12] + '_re_dedup_sort_realign.table'
args =  'java -Xmx4g -jar /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/GATK+Queue/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T BaseRecalibrator -R /Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta -I ' + realign_out + ' -o ' + recal_tbl
subprocess.Popen(args, shell=True).wait()


print 'done'
