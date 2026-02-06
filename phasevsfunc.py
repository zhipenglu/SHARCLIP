"""
phasevsfunc.py, zhipengluchina@gmail.com, 2025-07-19
checks HNRNPC arrays in the context of functional motifs, e.g. splice sites.
see phasevsOther.py for analysis of phasing for icSHAPE, eCLIP, motifs etc.
The data are still noisy, due to insufficient coverage in certain regions. 
Some starts and ends are not clearly defined from genomic studies. 

What's next?
1. incorporate density information
2. consider isoforms based on actual data.
3. include all samples. Some may have higher quality.
4. normalize intron vs. exon coverage. 

##### 0. Basic conclusions
A. Can we analyze HNRNPC occupancy?
B. There is a clear pattern, simply based on analysis of the peak cntrs
C. Annotations are not very accurate.
D. Need to consider cell type specific binding events. 
E. Do we need to consider the regularity score?

##### 1. Basic procedure
A. read annotations to get features: starts, ends, splices, Alu, etc. 
B. extract peaks in these regions.
C. average and plot the profiles

#bed12: chrom,start,end,name,score,strand,tstart,tend,rgb,bcount,bsizes,bstarts

##### 2. Example command
cd /Users/lu/Desktop/DG 
folder=/Users/lu/Documents/lulab/projects/hg38/refGene/
peaksbw=HepG2_HNRNPC_asso_gap1_sorted_shortcntrs.bw
python ~/Dropbox/_scripts/phasevsfunc.py ${folder}ACTB.bed $peaksbw c
python ~/Dropbox/_scripts/phasevsfunc.py ${folder}hg38refGene.bed $peaksbw c

"""



################################################################################
def getfeatures(bed):
    f = open(bed,'r') #each dict, value = gene name
    txstart,txend,tlstart,tlend = {},{},{},{} #transcription and translation
    ss5,ss3 = {},{} #splice sites, what about branchpoints?
    #key=(chrom,pos,strand), value=name
    for line in f:
        record = line.split()
        chrom,start,end,name,_,strand,tstart,tend,_,_,bsizes,bstarts = record
        start,end,tstart,tend = [int(i) for i in [start,end,tstart,tend]]
        #A. get splice sites first
        bcount = int(record[9])
        if bcount == 1: continue
        bsizes = [int(i) for i in bsizes.split(',')]
        bstarts = [int(i) for i in bstarts.split(',')]
        gstarts = [start+bstarts[i] for i in range(bcount)] #genomic starts/ends
        gends = [gstarts[i]+bsizes[i] for i in range(bcount)]
        for i in gends[:-1]: ss5[(chrom,i,strand)] = name
        for i in gstarts[1:]: ss3[(chrom,i,strand)] = name
        if strand == '-':
            gstarts,gends = gends,gstarts
            for i in gends[1:]: ss5[(chrom,i,strand)] = name
            for i in gstarts[:-1]: ss3[(chrom,i,strand)] = name
        #B. get transcription and translation starts/ends.
        if strand == '-': start,end,tstart,tend = end,start,tend,tstart
        txstart[(chrom,start,strand)] = name; txend[(chrom,end,strand)] = name
        if tstart == tend: continue #ncRNA
        tlstart[(chrom,tstart,strand)] = name; tlend[(chrom,tend,strand)] = name
    f.close()
    return txstart,txend,tlstart,tlend,ss5,ss3 
################################################################################





################################################################################
##### 2. overlay the HNRNPC peaks on genomic features, e.g. functional motifs
def featurecov(bw,features,r):
    #bw is the HNRNPC SHARCLIP peak centers
    #features, key=(chrom,pos,strand), value=name
    #r is radius from the feature point
    #covd: summary coverage in 200nt windows
    f = pyBigWig.open(bw,'r'); chroms = f.chroms()
    covd = {i:0 for i in range(-r,r)} #coverage dictionary
    for feature in features:
        chrom,pos,strand = feature
        if chrom not in chroms: continue
        wigs = f.intervals(chrom,max(0,pos-r),min(pos+r,chroms[chrom]))
        if not wigs: continue
        for w in wigs:
            k = w[0]-pos
            if strand == '-': k = -k
            if k in covd: covd[k] += 1
    return covd #key=pos, value=count
################################################################################





################################################################################
##### 3. plot the summary coverage
def plotcov(peaksbw,featurename,covd):
    figname = peaksbw[:-3]+"_around_"+featurename+".pdf"
    data = list(covd.values()); m = max(data); data = [i/m for i in data]
    plt.plot(data); plt.ylim(0,1); plt.savefig(figname); plt.close(); return
    
################################################################################





################################################################################
##### 4. process the input files
if __name__ == '__main__':
    import os; os.environ['OPENBLAS_NUM_THREADS'] = '1'
    import sys, pyBigWig, datetime
    import numpy as np
    from matplotlib import pyplot as plt

    if len(sys.argv) < 4:
        sys.exit("python phasevsfunc.py bed peaksbw regbw")
    bed,peaksbw,regbw = sys.argv[1:4]
    r = 100 #radius 100 on each side
    
    #A. make a set of dicts for features. key=(chrom,pos,strand), value=name
    print(str(datetime.datetime.today())+"\tMaking features dicts") 
    txstart,txend,tlstart,tlend,ss5,ss3 = getfeatures(bed)
    fd = {"txstart":txstart,"txend":txend,"tlstart":tlstart,"tlend":tlend,
          "ss5":ss5,"ss3":ss3}; #print("ss5", len(ss5))

    #B. extract data for these ranges.
    for d in fd: #d: feature name, e.g., txstart or ss5
        covd = featurecov(peaksbw,fd[d],r)
        plotcov(peaksbw,d,covd)     
        #for k in range(-r,r): print(covd[k])
    print(str(datetime.datetime.today())+"\tFinished analysis")
    
################################################################################

