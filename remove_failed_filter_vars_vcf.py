__author__ = 'mjohnpayne'


invcf = open("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/PM1_catchup/PM1_GATK_raw_out_process_snp_filt.vcf","r")
outvcf = open("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/PM1_catchup/PM1_GATK_snp_final.vcf","w")
for i in invcf:
    if i[0] == "#":
        outvcf.write(i)
    else:
        col = i.strip('\n').split('\t')
        if col[6] == "PASS":
            outvcf.write(i)
outvcf.close()