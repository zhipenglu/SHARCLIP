"""
phasevsOther.py, zhipengluchina@gmail.com, 2025-07-10
use xsite.py to identify crosslinking sites of RBPs.
export a summary over the arrays of peaks (only high regularity scores). 
then use plotphaseRBP.py to plot a heatmap of the coverage tracks.

#USE THESE ON SERVER:
conda activate env_full
import os; os.environ['OPENBLAS_NUM_THREADS'] = '1' 

#earlier steps, for the analysis of ribonucleosome organization:
1. positioning.py: short/long/all mid/density/positioning bws from gap1 bams
2. findpeaks.py: short/long/all peak cntrs.bw
2. peakperiodicity: summarize distributions and wave form. 
3. xsite.py: get crosslinking sites from eCLIP data
3. ntdensity.py: report kmer periodicity and positions. 
4. peakconnect.py: calculate connection among peaks using cntrs.bw and gap1.bam
4. phasediff.py: compare among sample. 
5. phaseregularity.py: loosepeaks.bw, loosepeaks_regularity.bw, hist.pdf/txt

Major comparisons:
1. RBP binding sites from eCLIP experiments
2. kmers, e.g. UUUU
3. splice sites.

Despite the lack of seq specificity, positioning was very consistent across
conditions, cell types, and likely represent the landscape where RBPs bind
Here are the questions. Where are they? How do they organize other RBPs
Does it explain why certain motifs are bound or not? Probably. 
icSHAPE data from Sun 2019 NSMB, quantity is low, so we will use SHARCLIP data
to guide ASO design.

For some of the low quality/quantity datasets, try merging replicates first.
For eCLIP for example, merge the 2 or 4 replicates of HepG2 and K562.
icSHAPE data, either from my own analysis, or the RASP database: 
folder=~/Documents/lulab/projects/shape/shape_subcell
hek_ch_vivo_minus_sorted.bw, hek_ch_vivo_plus_sorted.bw
hek_ch_vivo_imputed_neg.bw, hek_ch_vivo_imputed_pos.bw from RASP
covmin=0 for eCLIP due to lower quality and quantity. covmin=0.5 for icSHAPE.

Steps:
1. read peak centers bw as a list, only keeping ones with good regularity.
2. for each center position, extract SHAPE/CLIP in a window of several monomers.
3. if sufficient SHAPE/CLIP coverage, average data across the window
4. TODO: separate structures vs. monomer binding sites in future.
5. TODO: separate normal vs. repetitive regions later
6. TODO: separate exons vs. introns, may be differences
7. TODO: compare chromatin vs. later stages, e.g., cytoplasmic icSHAPE
8. TODO: separate asso vs. indep structrues ... 


comparing hnRNP phase with icSHAPE. 
################################################################################
##### Example using ACTB:
cd /Users/lu/Desktop/DG
genesbed=~/Documents/lulab/projects/hg38/refGene/hg38refGene_25chroms.bed
regbw=~/Desktop/DG/HepG2_HNRNPC_asso_gap1_sorted_short_loosepeaks_regularity.bw
cntrs=~/Desktop/DG/HepG2_HNRNPC_asso_gap1_sorted_shortcntrs_masked.bw
cntrs=~/Desktop/DG/HepG2_HNRNPC_asso_gap1_sorted_shortcntrs.bw
$folder/hek_ch_vivo_plus_sorted.bw, $folder/hek_ch_vivo_minus_sorted.bw
$folder/hek_cy_vivo_plus_sorted.bw,$folder/hek_cy_vivo_minus_sorted.bw
f=HepG2_HNRNPC_asso_gap1_sorted_shortcntrs.bw 
awk '$1=="chr7" && $2>5527148 && $3<5530601' f > ${f/bedgraph/_ACTB.bedgraph}
HepG2_HNRNPC_asso_gap1_sorted_shortcntrs_ACTB.bw
python ~/Dropbox/_scripts/phasevsOther.py $genesbed $cntrs $regbw $shapebw


comparing hnRNP phase with other RBPs. DONE for HepG2 and K562. 
################################################################################
##### 1. Comparing with HNRNPC eCLIP. HepG2_HNRNPC_0/1/2_xl.bw
genesbed=~/Documents/lulab/projects/hg38/refGene/hg38refGene_25chroms.bed
regbw=~/Desktop/DG/HepG2_HNRNPC_asso_gap1_sorted_short_loosepeaks_regularity.bw
cntrs=~/Desktop/DG/HepG2_HNRNPC_asso_gap1_sorted_shortpositioning_cntrs_masked.bw
cntrs=~/Desktop/DG/HepG2_HNRNPC_asso_gap1_sorted_shortpositioning_cntrs.bw
clipbw=~/Desktop/eCLIP_desktop/HepG2_HNRNPC_0_xl.bw
python ~/Dropbox/_scripts/phasevsOther.py $genesbed $cntrs $regbw $clipbw

genesbed=/project/zhipengl_72/labshare/hg38/hg38refGene_25chroms.bed
fo=/project/jianhuib_673/sharclip/sharclip_merged_gapped_reads/gap1/positioning/
regbw=${fo}HepG2_HNRNPC_asso_gap1_sorted_short_loosepeaks_regularity.bw
cntrs=${fo}HepG2_HNRNPC_asso_gap1_sorted_shortpositioning_cntrs.bw
python /project/zhipengl_72/labshare/scripts/phasevsOther.py \
$genesbed $cntrs $regbw K562_ZRANB2_2_xl.bw

K562_HNRNPC_asso_gap1_sorted_shortpositioning_cntrs.bw

"""


import os; os.environ['OPENBLAS_NUM_THREADS'] = '1' #to override a bug on CARC
import sys, pyBigWig, datetime
import numpy as np
from matplotlib import pyplot as plt
from struclib import readbed6






################################################################################
##### 1. all peaks are extracted now. Test subsets later. 
def getpeaks(peaksbw,regbw,genesdict,regmin):
    """
    extract peaks from bw, and assign strand information.
    only keep those with regularity score regmin >=0.3. 
    regmin is set to 0.3 based on visual inspection and distribution histogram.
    """
    peaksf = pyBigWig.open(peaksbw, 'r'); peaks = [] #[[chrom,pos,strand],...]
    regf = pyBigWig.open(regbw, 'r')
    chroms = peaksf.chroms()
    for gene in genesdict: #genesdict: {gene:[chrom,start,end,score,strand],...}
        chrom,start,end,score,strand = genesdict[gene]
        if chrom not in chroms: continue
        regwig = regf.intervals(chrom,start,end)
        if not regwig: continue
        for w1 in regwig:
            s,e,v = w1 #start,end,value
            if v < regmin: continue
            peakwig = peaksf.intervals(chrom,s,e)
            if not peakwig: continue
            for w2 in peakwig: peaks.append([chrom,w2[0],strand])
    peaksf.close(); regf.close()
    return peaks #[[chrom,pos,strand],...]
################################################################################





################################################################################
def rndpeaks(peaksbw,):
    """
    peaksbw contains prominences
    density file is useful for picking ...
    
    """
    return 
################################################################################





################################################################################
##### 2. count coverage for each bw file: 
def countall(bw,chroms24):
    f = pyBigWig.open(bw, 'r'); total = 0
    for chrom in chroms24:
        wig = f.intervals(chrom)
        if not wig: continue
        for w in wig: total += (w[1]-w[0])*w[2]
    return total
################################################################################






################################################################################
##### 3. peaks vs. SHAPE or RBP CLIP data bw files.
def vsother(peaks,bw,chroms24,radius,covmin,maxheight):
    #peaks: [[chrom,pos,strand],...]
    #SHAPE or RBP XL site file in bw format, both negative and positive files
    #chroms24: chr1-22,X,Y. chrM is ignored due to location.
    #chrM and MALAT2 are the extremely high regions in HNRNPC eCLIP. 
    #radius, around the peak center, default=100 ntss
    #covmin, min % of nts with a value in window (including 0). Depends on data
    #maxheight default 100, reduce impact of high peaks, likely artifacts.
    f = pyBigWig.open(bw,'r'); chroms = f.chroms(); covd = [] #list of lists
    for c in chroms24: chroms24[c] = chroms[c]
    for peak in peaks:
        chrom,pos,strand = peak
        if chrom not in chroms24: continue
        i1,i2 = max(pos-radius,0),min(pos+radius,chroms24[chrom])
        wig = f.intervals(chrom,i1,i2) #((start,end,value),...)
        if not wig or sum([w[1]-w[0] for w in wig])/(i2-i1) < covmin or \
           max([w[2] for w in wig]) >= maxheight: continue
        d = [0 for i in range(0,i2-i1)]
        for w in wig:
            for j in range(w[0],w[1]): d[j-pos+radius] = w[2]
        if strand == '-': d = d[::-1]
        covd.append(d)
    f.close()
    return covd #list of lists, to be averaged later. 
################################################################################





    
################################################################################
##### 5. process the input files. 
if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("python phasevsOther.py genesbed peaksbw regbw bwfiles")
    radius,covmin = 100, 0 ##### 0 for eCLIP, 0.5 for SHAPE.
    chroms24 = {"chr"+i:0 for i in[str(j)for j in list(range(1,23))]+["X","Y"]}
    maxheight = 100 #may need to change for datasets other than HNRNPC. 
    genes,peaksbw,regbw,bw = sys.argv[1:5]; regmin = 0.4

    #A. Get the genes dictionary 
    print(str(datetime.datetime.today())+"\tGetting genes dictionary") 
    genesdict = readbed6(genes)

    #B. Only use good peaks, with high regularity scores
    print(str(datetime.datetime.today())+"\tGetting SHARCLIP peaks") 
    peaks = getpeaks(peaksbw,regbw,genesdict,regmin); #print(len(peaks))
    #3354648 peaks total, 455335 with regmin=0.3, 146220 with regmin=0.4

    #C. Compute average coverage over the windows 
    print(str(datetime.datetime.today())+"\tComputing cov for sample 1")
    covd1 = vsother(peaks,bw,chroms24,radius,covmin,maxheight)
    avg1 = np.mean(np.array(covd1),axis=0)
    outf1 = open(peaksbw[:-3]+'_'+bw.split('/')[-1][:-3]+".txt", 'w')
    outf1.write(str(countall(bw,chroms24))+'\n')
    outf1.write(''.join([str(i)+'\n' for i in avg1])); outf1.close()
    print(str(datetime.datetime.today())+"\tFinished analysis") 

################################################################################






