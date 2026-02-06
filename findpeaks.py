"""
findpeaks.py, zhipengluchina@gmail.com, 2025-07-04
input: positioning.bw, output: positioning_cntrs.bw, value = peak prominence
A wrapper for function peakfinder() in struclib.py. Default parameters:
{"height":0.2,"distance":30,"prominence":0.05,"width":(5,40),"wlen":100}

Major steps:
1. positioning.py->12 files, short/long/all/shortcx mid/density/positioning bw
2. Run this script using a slurm for array jobs: findpeaks.slurm

for occasional errors in bedgraph files,
use uniq -u to remove the duplicate lines. 

##### Example command for SHARCLIP and eCLIP data:
script=~/Dropbox/_scripts/findpeaks.py
bg2bw=bedGraphToBigWig
genome=~/Documents/lulab/projects/hg38/hg38.chrom456.sizes
bw=HepG2_HNRNPC_asso_gap1_PICALM_shortpositioning.bw
bw=HepG2_HNRNPC_0_positioning.bw
python $script $bw $bg2bw $genome

##### on CARC
script=/project/zhipengl_72/labshare/scripts/findpeaks.py
bg2bw=/project/zhipengl_72/labshare/kentutils/bedGraphToBigWig
genome=/project/zhipengl_72/labshare/hg38/hg38_SHARCLIP_sizes.txt
bw=HepG2_HNRNPC_asso_gap1_sorted_allpositioning.bw

for f in *allpositioning.bw; \
do (salloc python $script $f $bg2bw $genome &); done
ls -hl *all*cntrs*

python $script $bw $bg2bw $genome
"""

import os; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys, os, pyBigWig, datetime
from struclib import peakfinder

if __name__ == "__main__":    
    if len(sys.argv)<4: sys.exit("python findpeaks.py inbw bg2bw genome")
    bw = sys.argv[1]; chunk = 1000 #default, does not affect results.
    cntrsbg = bw[:-3]+"_cntrs.bedgraph"
    cntrsbw = bw[:-3]+"_cntrs.bw"
    uniq = 5 #only keep peaks if >5 uniq values in 20nt window.
    
    #1. export a bedgraph file for peak centers.
    print(str(datetime.datetime.today())+"\tFinding peaks")
    peakfinder(bw,chunk,cntrsbg,[],uniq) #default parameters set in struclib.py
    
    #2. convert the bedgraph to bw. 
    print(str(datetime.datetime.today())+"\tConverting bedgraph to bw")
    cmd = sys.argv[2] #absolute path of bg2bw
    genome = sys.argv[3] #genome size file
    line = "${} ${} ${} ${}".format(cmd,cntrsbg,genome,cntrsbw)
    print(line)
    os.system(line); #os.system("rm {}".format(cntrsbg))
    print(str(datetime.datetime.today())+"\tFinished analysis")

