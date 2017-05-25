import subprocess
import shlex
import sys
import re
import time
import glob
import os



t0= time.time()

def time_out():
    t= time.time()
    diff = t - t0
    timer = float(diff)
    timer = timer/60
    print '%1.2f'%timer + ' minutes'


inp = sys.argv[1]

file_lis = glob.glob(inp + '/*.fastq')
file_lis.sort()

for f in file_lis:
    print f.split('/')[-1].strip('.fastq')
a = 0
b = 1

while b < len(file_lis):
    f1 = file_lis[a]
    f2 = file_lis[b]

    print '\nRunning ' + f1.split('/')[-1].strip('_trim_1.fastq')
    #bwa mem -M -t 7 ref.fa read1.fq read2.fq > aln.sam

    sam_path = f1.strip('_trim_1.fastq')+".sam"

    sampe_args = '/usr/local/bin/bwa mem -M -t 7 /Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13.fasta ' + f1 + ' ' + f2 + ' ' + ' > ' + sam_path

    sampe_out = subprocess.Popen(sampe_args, shell=True).wait()

    print '\n\n'

    time_out()

    bam_path = sam_path.strip('sam') + 'bam'

    to_bam_args = '/usr/local/bin/samtools view -S -b -o ' + bam_path + ' ' + sam_path

    print '\n\nConverting sam to bam file\n'

    bam_out = subprocess.Popen(to_bam_args, shell=True).wait()

    print '\n\n'

    time_out()

    sort_bam_path = bam_path.strip('.bam') + '_sort'

    sort_bam_args = '/usr/local/bin/samtools sort ' + bam_path + ' ' + sort_bam_path

    print '\n\nSorting bam file\n'

    bam_sort = subprocess.Popen(sort_bam_args, shell=True).wait()

    index_args = '/usr/local/bin/samtools index ' + sort_bam_path

    subprocess.Popen(index_args, shell=True).wait()

    print '\n\n'

    os.remove(sam_path)
    os.remove(bam_path)

    time_out()
    a+=2
    b+=2


