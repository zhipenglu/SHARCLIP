"""
bam2heat.py, Zhipengluchina@gmail.com. 2025-07-01. 
updated 2025-09-03: random subsample instead of uniform scale for normalization
updated 2025-10-10: plotting composite samples. 
Either plots pairs of individual samples, or composite samples, e.g.
HNRNPCassomul vs. HNRNPCindepmul. The composite samples increases robustness
for regions with low coverage, but does not allow comparisons across samples. 

see bam2DGorg.py for additional analysis of high level organizations
see processDG.py for plotting DGs as arcs with coverage scaling
see tgassembly.py for details on the usage of intervaltree
major samples: 7 cell lines, 5 RBPs, asso and indep, and other
Also check other samples from published studies for comparison and validation

##### all purposes
1. normalize between samples. 
2. plot heatmap for coverage along and off the diagonal, for one sample. 
3. plot heatmap for single samples or in pairs, or difference between the two.


################################################################################
##### Examples plotted so far: 
ATAD3B: s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
python $script $bam $s1 $s2 ATAD3B chr1 1471500 32000 20 all #1477000 11000 5

EXOSC6: s1=HepG2HNRNPCassomul; s2=HepG2HNRNPCindepmul
chr16:70246778-70251940 #EXOSC6 chr16 70246778 5200 5

HNRNPLL: s1=HepG2HNRNPCassomul; s2=HepG2HNRNPCindepmul
38560000 44000 100 #38564000 5000 5 #38565300 1800 5 #38594500 2200 5

RBFOX2: s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
RBFOX2 chr22 35738500 30000 20 #35751500 9000 5

RPS10. chr6:34417000-34426500 #RPS10 5 9500 34417000
s1=HepG2HNRNPCassomul s2=HepG2HNRNPCindepmul

NUMB. classic example of SE; s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
NUMB 73275000 183800 100 #73313000 7000 20 #73310000 14000 50

ZNF276; s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
python $script $s1 $s2 ZNF276 20 20000 89721000 # 5 3000 89732000

DNM1, clear switching in development
chr9:128203000-128256000. #DNM1 128222000 13000 100
s1=HepG2HNRNPCassomul; s2=HepG2HNRNPCindepmul
s1=HepG2HNRNPCassomul; s2=HepG2HNRNPCassomul
s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
python $script $bam $s1 $s2 DNM1 chr9 128228000 1000 5 all

s1=astroHNRNPCassomul; s2=neuronHNRNPCassomul
s1=astroHNRNPCindepmul; s2=neuronHNRNPCindepmul
python $script $bam $s1 $s2 DNM1 chr9 128222000 13000 30 all

DNM2. classic example of developmental switching. 
s1=HepG2HNRNPCassomul; s2=HepG2HNRNPCindepmul
python $script $s1 $s2 DNM2 chr 10717900 114300 100 #10794000 19000 100

ZSCAN25. chr7:99616000-99740000
s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
ZSCAN25 99616000 124000 300 #99616500 16200 20 #99622000 5000 5

ENAH, example from Lovci 2013.
s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
s1=K562HNRNPCassomul; s2=K562RBFOX2assomul
s1=K562RBFOX2assomul; s2=K562RBFOX2indepmul
python $script $bam $s1 $s2 ENAH chr1 225500500 8000 50 all #225503000 2100 20
python $script $bam $s1 $s2 ENAH chr1 225503000 2100 5 all

KIF21A. chr12:39326000-39331000
s1=iPSCHNRNPCassomul; s2=neuronHNRNPCassomul
s1=iPSCHNRNPCindepmul; s2=neuronHNRNPCindepmul
s1=HNRNPCassomul; s2=HNRNPCindepmul
python $script $bam $s1 $s2 KIF21A chr12 39326000 5000 10 all

KIF2A. Zhang 2016. A5SS: chr5:62,352,400-62,353,585; SE: chr5:62366400-62373900
s1=neuronHNRNPCassomul; s2=astroHNRNPCassomul
s1=neuronHNRNPCindepmul; s2=astroHNRNPCindepmul
python $script $bam $s1 $s2 KIF2A chr5 62366400 7500 30 all

KIF1B. Clear example of RBFOX2 associated loop
s1=HepG2HNRNPCassomul; s2=HepG2HNRNPCindepmul
s1=HepG2HNRNPCassomul; s2=HepG2RBFOX2assomul
s1=HepG2RBFOX2assomul; s2=HepG2RBFOX2indepmul
s1=neuronHNRNPCassomul; s2=astroHNRNPCassomul
python $script $bam $s1 $s2 KIF1B chr1 10278500 4000 50 all #10272000 11000 30

PKM, classic example of MXE;
s1=neuronHNRNPCassomul; s2=iPSCHNRNPCassomul
s1=neuronHNRNPCindepmul; s2=iPSCHNRNPCindepmul
python $script $bam $s1 $s2 PKM chr15 72200000 3500 10 all
python $script $bam $s1 $s2 PKM chr15 72200000 7000 10 all


MTSS1, seems to have strong structure differences among cell types
s1=neuPCHNRNPCassomul; s2=neuronHNRNPCassomul
s1=astroHNRNPCassomul; s2=neuronHNRNPCassomul
python $script $bam $s1 $s2 MTSS1 chr8 124553000 10000 100 all

ADD1.strong SE. Reported by Zhang Xiaochang 2016. Qi 2021 Same lab.
One of the best examples so far. 
s1=astroHNRNPCassomul; s2=neuronHNRNPCassomul #astro is 3x deeper than neuron. 
s1=periHNRNPCassomul; s2=neuronHNRNPCassomul  #peri is only 10% deeper. 
s1=neuPCHNRNPCindepmul; s2=neuronHNRNPCindepmul 
python $script $bam $s1 $s2 ADD1 chr4 2914000 16000 50 all #2926000 2200 10
python $script $bam $s1 $s2 ADD1 chr4 2926000 1200 5

FLNA, strong SE. Two regions; s1=iPSCHNRNPCassomul; s2=astroHNRNPCassomul 
s1=iPSCHNRNPCindepmul; s2=astroHNRNPCindepmul 
python $script $bam $s1 $s2 FLNA chrX 154354800 3000 10 all #154355130 2300 10

DLG1, Zhang 2016. Strong switching across 7 cell types. 
s1=periHNRNPCassomul; s2=neuronHNRNPCassomul 
s1=periHNRNPCindepmul; s2=neuronHNRNPCindepmul 
python $script $bam $s1 $s2 DLG1 chr3 197069000 8000 30 all #197074300 2400 10 

SCRIB, Zhang 2016. Good series of SE. gap1 cov 20-50, not very high. 
s1=astroHNRNPCassomul; s2=neuronHNRNPCassomul 
python $script $bam $s1 $s2 SCRIB chr8 143791440 430 5 trim

TPM2, strong switching in development. data not enough for some samples
s1=astroHNRNPCassomul; s2=neuronHNRNPCassomul
s1=astroHNRNPCindepmul; s2=neuronHNRNPCindepmul
s1=iPSCHNRNPCassomul; s2=astroHNRNPCassomul

python $script $bam $s1 $s2 TPM2 chr9 35684400 1000 5 all
Analysis reveals a strong duplex covering one of the two MXEs, switching is not
clear yet. Needs deeper analysis. 


EEF1D long exon; s1=K562HNRNPCassomul; s2=neuronHNRNPCassomul
python $script $bam $s1 $s2 EEF1D chr8 143586500 11000 30 all #143587500 2000 5

ACTN4, MXE, strong differences, clear switch mechanism. 
s1=astroHNRNPCassomul; s2=neuronHNRNPCassomul
s1=astroHNRNPCindepmul; s2=neuronHNRNPCindepmul
python $script $bam $s1 $s2 ACTN4 chr19 38710000 1500 5 all #38709300 2400 5

SNX21; s1=neuPCHNRNPCassomul; s2=neuronHNRNPCassomul
python $script $bam $s1 $s2 SNX21 chr20 45834500 6500 30 all

CTTN; s1=iPSCHNRNPCassomul; s2=periHNRNPCassomul #also indep
s1=iPSCHNRNPCindepmul; s2=K562HNRNPCindepmul
s1=iPSCHNRNPCindepmul; s2=astroHNRNPCindepmul
python $script $bam $s1 $s2 CTTN chr11 70420400 2600 10 all

SCARB1. two groups of high and low inclusion among 5 cell types. Good gap1 cov. 
s1=iPSCHNRNPCindepmul; s2=neuPCHNRNPCindepmul
s1=iPSCHNRNPCindepmul; s2=HepG2HNRNPCindepmul
python $script $bam $s1 $s2 SCARB1 chr12 124777600 9000 10 all

ATE1: at most 3 fold switching, including HepG2 and K562. 
s1=K562HNRNPCassomul; s2=astroHNRNPCassomul #also indep
s1=K562HNRNPCindepmul; s2=astroHNRNPCindepmul #also indep
s1=K562HNRNPCassomul; s2=neuPCHNRNPCassomul #also indep
s1=K562HNRNPCindepmul; s2=neuPCHNRNPCindepmul #also indep
s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul #also indep
s1=iPSCHNRNPCassomul; s2=neuPCHNRNPCassomul #also indep
s1=iPSCHNRNPCindepmul; s2=neuPCHNRNPCindepmul #also indep

python $script $bam $s1 $s2 ATE1 chr10 121898800 1400 5 all

python $script $bam $s1 $s2 ATE1 chr10 121870000 33000 100 all


CASK. Experimentally validated by Margasyuk et al. 2023
s1=neuronHNRNPCassomul; s2=neuPCHNRNPCassomul
s1=neuronHNRNPCindepmul; s2=neuPCHNRNPCindepmul
s1=neuronHNRNPCindepmul; s2=HepG2HNRNPCindepmul
python $script $bam $s1 $s2 CASK chrX 41553500 6500 10 all

PHF21A, very strong MXE, good signal in all 7 lines. gap1 cov 130-320. 
s1=neuronHNRNPCassomul; s2=periHNRNPCassomul
python $script $bam $s1 $s2 PHF21A chr11 45938000 12000 50 all #45945500 3500 10

BRD2/BRD3. no difference among cell lines, but validated by ASO. Petrova 2024
s1=iPSCHNRNPCassomul; s2=astroHNRNPCassomul
s1=iPSCHNRNPCindepmul; s2=astroHNRNPCindepmul
python $script $bam $s1 $s2 BRD2 chr6 32974750 650 5 all

AGRN, strong inter-intron contacts without skipping.
s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul
python $script $bam $s1 $s2 AGRN chr1 1051700 3000 10 all

TNFRSF25, no clear AS; s1=neuronHNRNPCassomul; s2=astroHNRNPCassomul #also indep
python $script $bam $s1 $s2 TNFRSF25 chr1 6464600 900 5 all

CLSTN1, another with strong inter-intron duplexes but no skipping
there are even structures that block splice sites. these are normal, but why?
s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul
python $script $bam $s1 $s2 CLSTN1 chr1 9730700 1100 5 all

Other good examples from Zhang 2016 Cell. Not plotted yet. 
CLIP1, chr12:122,269,296-122,283,228, neuron specific, conserved. Took picture
CLIP2, chr7:74,371,826-74,375,989, Took picture
CLASP1, chr2:121406035-121412364, high cov. several SE and MXE. Took picture
CLASP2, At least 5 SEs. Strong diff. Took picture. 
MAST2, chr1:45,994,726-46,007,738, SE and RI. broad conservation. Took picture. 
MARK4, neuron specific, good coverage. took picture. 
DYNC1I2, neuron specific, broad conservation, high coverage. Took picture. 
DCTN1, 3 neuron selective exons, broad conservation. Took picture
EPB41L1, cov good, spread in 7 lines, two SE. Took picture. 
EPB41L3, two groups of AS, neuron selective. Took picture. 
SYNE2, chr14:64,213,958-64,216,754, good coverage, consered, 
DCLK1, HOMER1, GPHN, no clear AS

ZKSCAN1, circular RNA; s1=astroHNRNPCassomul; s2=astroHNRNPCindepmul
python $script $bam $s1 $s2 ZKSCAN1 chr7 100015000 15000 100 all
python $script $bam $s1 $s2 ZKSCAN1 chr7 100022800 2400 10 all

HIPK3, circular RNA; s1=astroHNRNPCassomul; s2=astroHNRNPCindepmul
python $script $bam $s1 $s2 HIPK3 chr11 33257000 72000 200 all #33285500 3000 10

EPHB4, circular RNA; s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul
python $script $bam $s1 $s2 EPHB4 chr7 100807300 6500 50 all

FOXP1, critical regulator of stem cell TF network. Gabut ... Blencowe 2011.
Not easy to identify the switching mechanisms. Too many DGs, many associated
with splicing, but not good coverage. 
s1=iPSCHNRNPCassomul; s2=neuPCHNRNPCassomul 
python $script $bam $s1 $s2 FOXP1 chr3 70965000 10000 30 all
python $script $bam $s1 $s2 FOXP1 chr3 70968000 10000 30 all

s1=iPSCHNRNPCindepmul; s2=astroHNRNPCindepmul 
python $script $bam $s1 $s2 FOXP1 chr3 70970000 7200 5 all
s1=iPSCHNRNPCassomul; s2=periHNRNPCassomul 
s1=iPSCHNRNPCindepmul; s2=periHNRNPCindepmul 
python $script $bam $s1 $s2 FOXP1 chr3 70970000 7200 5 all
python $script $bam $s1 $s2 FOXP1 chr3 70971500 1500 5 all

s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul 
python $script $bam $s1 $s2 SLC25A3 chr12 98593500 5000 20 all

Google search suggests QKI as another key regulator of development as a result
of its own splicing. Our 7 cell lines have sufficient data to call the AS and
also enough gap1 reads to analyze structures.

PTBP, CELF, RBFOX and MBNL family of RBPs are also heavily alternatively spliced
to control downstream splicing programs. Clear SE in MBNL1/2/3, Clear RI in
CELF1


TCF3, an essential regulator of stem cells, also has an RI next to it and A5SS.
One of the best examples, in terms of function and structural mechanisms. 
See three studies, Yamazaki 2018, 2019 and 2020. Earlier work by ... 
s1=iPSCHNRNPCassomul; s2=periHNRNPCassomul #high coverage, a bit slow. 
python $script $bam $s1 $s2 TCF3 chr19 1609000 11000 30 all

s1=iPSCHNRNPCassomul; s2=periHNRNPCassomul 
s1=iPSCHNRNPCindepmul; s2=periHNRNPCindepmul 
python $script $bam $s1 $s2 TCF3 chr19 1612000 3600 10 all
python $script $bam $s1 $s2 TCF3 chr19 1611500 4500 10 all

TCF7L2
s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul 
s1=HepG2HNRNPCindepmul; s2=periHNRNPCindepmul
s1=neuronHNRNPCindepmul; s2=periHNRNPCindepmul
python $script $bam $s1 $s2 TCF7L2 chr10 113158000 10000 30 all

s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul 
s1=HNRNPCassomul; s2=HNRNPCindepmul 
python $script $bam $s1 $s2 SMN2 chr5 70049600 28000 100 all #c 100 50
python $script $bam $s1 $s2 SMN2 chr5 70070600 7000 20 all #c 20 10
python $script $bam $s1 $s2 SMN2 chr5 70076300 1300 5 all #c 20 2
python $script $bam $s1 $s2 SMN2 chr5 70076500 600 5 all #c 5 1

##### should not use all. Use other options to increase resolution. 

s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul
s1=HepG2RBFOX2assomul; s2=HepG2RBFOX2assomul
python $script $bam $s1 $s2 MYH10 chr17 8569000 9000 20 all
python $script $bam $s1 $s2 MYH10 chr17 8574400 3000 10 all

s1=astroHNRNPCindepmul; s2=periHNRNPCindepmul 
s1=astroHNRNPCassomul; s2=astroHNRNPCindepmul 
python $script $bam $s1 $s2 COL5A1 chr9 134825000 8000 30 all #c 20 2
python $script $bam $s1 $s2 COL5A1 chr9 134829900 500 5 all #c 20 2

s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul 
python $script $bam $s1 $s2 SLC39A14 chr8 22408000 7000 20 all
python $script $bam $s1 $s2 SLC7A2 chr8 17554000 2000 5 all

s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul 
s1=astroHNRNPCassomul; s2=astroHNRNPCindepmul 
python $script $bam $s1 $s2 ANKRD36C chr2 95876000 66000 100 all #c 200 100
python $script $bam $s1 $s2 ANKRD36C chr2 95895200 8000 20 all #50 20

IDE insulin degrading enzyme:
python $script $bam $s1 $s2 IDE chr10 92475500 8000 20 all # 50 20
python $script $bam $s1 $s2 IDE chr10 92478200 2200 5 all # 10 5

PTBP2, a strong SE and one of the most conserved introns.
s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul 
python $script $bam $s1 $s2 IDE chr1 96804500 3000 10 all # 50 20
very clear long range structures covering the SE, should be testable using ASOs

DPF2, Nazim 2024, Black lab, regulates stem cell and neuro diff.
s1=iPSCHNRNPCassomul; s2=neuronHNRNPCassomul 
s1=iPSCHNRNPCindepmul; s2=neuronHNRNPCindepmul 
s1=HNRNPCassomul; s2=HNRNPCindepmul

s1=-RIC; s2=-RIC
bam=sharclip_ric_cric_eclip_karr_sharc_paris_all_gap1_DPF2.bam
python $script $bam $s1 $s2 DPF2 chr11 65333800 21000 50 all # 50 20
python $script $bam $s1 $s2 DPF2 chr11 65343700 2700 10 all # 50 20

python $script $bam $s1 $s2 ENO1 chr1 8868000 4000 10 all 

s1=iPSCHNRNPCindepmul; s2=neuronHNRNPCindepmul;
python $script $bam $s1 $s2 FLOT2 chr17 28883000 6000 20 all 
python $script $bam $s1 $s2 FLOT2 chr17 28883800 2400 10 all 


s1=astroHNRNPCindepmul; s2=neuronHNRNPCindepmul;
python $script $bam $s1 $s2 RAB6A chr11 73716000 5000 20 all 

UTY. combination of multiple-A3SS and SE
s1=iPSCHNRNPCassomul; s2=iPSCHNRNPCindepmul 
python $script $bam $s1 $s2 UTY chrY 13357800 2200 5 all 

##### common variables:
cd ~/Documents/lulab/projects/sharclip; mamba activate crssant;
script=~/Dropbox/_scripts/bam2heat.py
bam=sharclip_ric_cric_eclip_karr_sharc_paris_all_gap1_mRNAs_Birch_\
T30_sorted_SHARCLIP.bam; #does not have to match other files BIRCH parameters.  
python $script $bam $s1 $s2 name chrom start length binsize all/trim/median

python $script $bam $s1 $s2 SMN2 chr5 70076500 600 5 all #c 5 1

"""


    

################################################################################
##### 1. parameters to be set in this script.
def setparams():
    #bothc and diffc are set empirically each time to draw legible heatmaps
    #use bigger numbers for zoom out views, and smaller for zoom in.
    #typically diffc is set as half or a quarter of bothc
    #For example in DNM1, full length: bothc=1000, 13kb: 200, 1kb: 20, 200nt: 2
    
    bothc = 50 #norm color for both triangle heatmaps, higher more saturation
    diffc = 25 #norm differential heatmap color. Higher means more saturation    radius = 7 #no need to change, radius around the midpoints
    matrixnorm = 0.5 #0.5 between RBPs, e.g. RBFOX2 vs HNRNPC; 0.9 asso vs indep
    iden = 0 #indicator for plotting two identical samples.
    spanmax = 1E9 #plot shortest span alignments. Use 1E9 to disable
    return bothc,diffc,matrixnorm,iden,spanmax
################################################################################





################################################################################
##### 2. export to heatmap using a scatter plot method, minimal pdf size
def heatmapscatter(start,rnasize,binsize,matrix1,matrix2,heatmap,bothc,iden):
    # A. plot the heatmap. Currently only using this option
    fig, ax = plt.subplots(); plt.gca().set_aspect('equal')
    s,b = start,binsize; ncol = int(rnasize/b+1); r = list(range(ncol)); d = {}
    for i,j in product(r,r):
        if i<=j:
            if matrix1[i][j]: d[(j,i)] = matrix1[i][j]
            if matrix2[i][j]: d[(i,j)] = matrix2[i][j]
    maxv = max(list(d.values()))
    for k in d.keys():
        c = min(1,d[k]/maxv*bothc) #scaled intensity.
        if not iden: c = (1,0,0,c) if k[0]>k[1] else (0,0,1,c)
        else: c = str(1-c) #diagonal coverage
        #c = (1-c/2,1-c,1-c/2) if iden else str(1-c), #dark purple
        ax.add_patch(patches.Rectangle((k[0],k[1]),1,1,facecolor=c))
    plt.xlim(0,ncol); plt.ylim(0,ncol); plt.savefig(heatmap); plt.close()
    print("Highest coverage in this heatmap:",maxv)
    return
################################################################################





################################################################################
##### 3. export to differential heatmap using a scatter plot method
def heatmapdiff(rnasize,binsize,matrixd,heatmap,diffc):
    fig, ax = plt.subplots(); plt.gca().set_aspect('equal')
    ncol = int(rnasize/binsize+1); r = range(ncol)
    d = dict(((j,i), matrixd[i][j]) for i in r for j in r if matrixd[i][j])
    minv,maxv = matrixd.min(),matrixd.max() #normalize against whole transcript
    maxb = min(abs(minv),abs(maxv)) #symmetric
    for k in d.keys():
        a = min(abs(d[k])/maxb*diffc,1); c = (1,0,0,a) if d[k]>0 else (0,0,1,a)
        ax.add_patch(patches.Rectangle((k[0],k[1]),1,1,facecolor=c))
        """
        #plot triangle for each contact. Not very useful. 
        if d[k]<0: continue #only plot the positive values. 
        c = (c[:3],c[3]/50); co1,co2 = k[0]+0.5,k[1]+0.5
        ax.add_line(lines.Line2D([co1,co1],[co1,co2],color=c,lw=1))
        ax.add_line(lines.Line2D([co1,co2],[co2,co2],color=c,lw=1))
        """
    plt.xlim(0,ncol); plt.ylim(0,ncol); plt.savefig(heatmap); plt.close();
    return
################################################################################





################################################################################
##### 4. process data for plotting.
if __name__ == "__main__":
    import sys, datetime
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches; from matplotlib import lines
    from itertools import product;
    from struclib import bam2matrix,bam2matrixSub

    if len(sys.argv) < 10: 
        sys.exit("bam2heat bam s1 s2 gene chrom start rnasize binsize option")
    sam,s1,s2,gene,chrom = sys.argv[1:6] #sample1 and sample2
    start,rnasize,binsize = int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8])
    radius = 7 #only used for the trim option
    option = sys.argv[9] #'all', 'trim' or 'median'
    list1 = (sam[:-4],s1,s2,start,rnasize,binsize)
    name1 = '{}_{}_{}_{}_{}_{}nt_heatmap_'.format(*list1)
    heatmapb,heatmapd = name1+'both.pdf',name1+'diff.pdf'
    heatmapcov1,heatmapcov2 = name1+'cov1.pdf',name1+'cov2.pdf' #on diagonal
    bothc,diffc,matrixnorm,iden,spanmax = setparams()
    
    ##### A. make cov matrices from bam file and normalize
    print(str(datetime.datetime.today())+"\tGenerating coverage matrices")
    #a = (sam,s1,s2,radius,chrom,start,rnasize,binsize,spanmax,matrixnorm,'all')
    #M0,M1,Ms = bam2matrix(*a) #use 90th percentile for normalization of matrix
    b = (sam,s1,s2,radius,chrom,start,rnasize,binsize,spanmax,'all')
    M0,M1,Ms = bam2matrixSub(*b) #use subsampling for normalization of reads

    ##### B. 
    ss = ['HNRNPCassomul','HNRNPCindepmul']; n = int(rnasize/binsize+1)
    c1,m1,ms1 = np.zeros((n,n)),np.zeros((n,n)),np.zeros((n,n));
    c2,m2,ms2 = np.zeros((n,n)),np.zeros((n,n)),np.zeros((n,n));
    md = m1-m2; msd = ms1-ms2
    if s1 in ss and s2 in ss: #s1/s2 are composite, HNRNPC[asso/indep]mul
        for k in M0: #k is the same across M0, M1 and Ms
            if s1 in k: c1 += M0[k]; m1 += M1[k]; ms1 += Ms[k]
            elif s2 in k: c2 += M0[k]; m2 += M1[k]; ms2 += Ms[k]
        md = m1-m2; msd = ms1-ms2
    else: #s1/s2 are individual samples, e.g. astroHNRNPC[asso/indep]mul
        c1,c2 = M0[s1],M0[s2] #coverage along diagonal
        m1,m2 = M1[s1],M1[s2]; md = m1-m2 #all contacts and diff
        ms1,ms2 = Ms[s1],Ms[s2]; msd = ms1-ms2 #short range contacts and diff
        #here spanmax was set to 1E9, so the short range ones are the full set.
        

    ##### C. plot heatmaps. 
    print(str(datetime.datetime.today())+"\tPlotting heatmaps and graphs")
    heatmapscatter(start,rnasize,binsize,ms1,ms2,heatmapb,bothc,iden)
    heatmapdiff(rnasize,binsize,md,heatmapd,diffc) 
    #bothc=2; iden=1; #plotting two identical samples, only for 
    #heatmapscatter(rnasize,binsize,c1,c1,heatmapcov1)
    #heatmapscatter(rnasize,binsize,c2,c2,heatmapcov2)
    print(str(datetime.datetime.today())+"\tFinished")
################################################################################





