__author__ = 'mjohnpayne'

import subprocess
import shlex
import sys
import re
import time
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import os
from datetime import datetime
from time import sleep as sl

my_env = os.environ
my_env["PATH"] = "/opt/local/bin:/opt/local/sbin:" + my_env["PATH"]

fasta = sys.argv[1]

var = sys.argv[2]

sss = '/Volumes/MP_HD/bys_snpeff_catchup/provean_snp_inputs/' + fasta.split('/')[-1].replace('.fasta','.sss')

out = '/Volumes/MP_HD/bys_snpeff_catchup/provean_snp_inputs/' + fasta.split('/')[-1].replace('.fasta','.out')

protean_args = '/usr/local/bin/provean.sh -q ' + fasta + ' -v ' + var + ' --save_supporting_set ' + sss
outp = open(out,'w')
print fasta
subprocess.Popen(protean_args, shell=True,stdout=outp,env=my_env).wait()
outp.close()