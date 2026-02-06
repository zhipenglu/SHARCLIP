"""
peakperiodicity.py, zhipengluchina@gmail.com, 2025-07-04.
input: peak center bw files. output: interpeak and peakperiod pdf and txt files
1. plots inter-peak distance distribution, for all peaks and tight clusters
2. plots periodicity of peaks, for all peaks and tight clusters. 

Note: Process all cntrs.bw files. SHARCLIP data for other proteins also have
this periodicity, likely for two reasons, enrichment on ribonucleosomes, or
their own periodic binding. Also check ribosome proteins on rRNAs and mature
mRNAs. The 39nt periodicity is to be present for other RBPs as well, in both
SHARCLIP and eCLIP, confirming the ribonucleosome as a fundamental unit of
assembly. The SHARCLIP data may be complicated due to extensive
protein-protein crosslinking. The periodicity from eCLIP strongly suggests
periodic patterns.

previous steps before this step:
1. positioning.py, extract midpoints, density and positioning from gap1 bam
2. findpeaks.py, get peak centers from the short/long/all positioning.bw

##### Major considerations:
1. Peak centers are used since midpoints, density, positioning are hard to align
2. Start with 400nt windows. Peaks decay within this range.
3. Tight clusters show better periodicity than all peaks, as expected.

##### small test 
python ~/Dropbox/_scripts/peakperiodicity.py PICALM_cntrs.bw tight

##### SHARCLIP data, HepG2, HNRNPC, whole genome, without masking
python ~/Dropbox/_scripts/peakperiodicity.py \
HepG2_HNRNPC_asso_gap1_sorted_shortcntrs.bw tight # in 10 sec. 

##### SHARCLIP data, HepG2, HNRNPC, whole genome, with masking,
bedtools subtract -a HepG2_HNRNPC_asso_gap1_sorted_shortcntrs.bedgraph -b \
/Users/lu/igv/genomes/hg38/hg38_mask_sorted.bed > \
HepG2_HNRNPC_asso_gap1_sorted_shortcntrs_masked.bedgraph #reduced 132 to 83M
python ~/Dropbox/_scripts/peakperiodicity.py \
HepG2_HNRNPC_asso_gap1_sorted_shortcntrs_masked.bw tight

##### For eCLIP data, whole genome
cd /Users/lu/Desktop/eCLIP_desktop
python ~/Dropbox/_scripts/peakperiodicity.py HepG2_HNRNPC_1_cntrs.bw tight

##### eCLIP data, after masking: 
bedtools subtract -a HepG2_HNRNPC_1_cntrs.bedgraph -b \
/Users/lu/igv/genomes/hg38/hg38_mask_sorted.bed > \
HepG2_HNRNPC_1_cntrs_masked.bedgraph
python ~/Dropbox/_scripts/peakperiodicity.py \
HepG2_HNRNPC_1_cntrs_masked.bw tight
"""





################################################################################
##### 1. calculate the decaying wave distribution. 
def peakperiod(bw,nsubset,spanmax,tight):
    """
    bw: peak centers file, value=peak prominence, but not used here. 
    nsubset: num of consecutive peaks to pick in each window. 
    spanmax: region length to plot the histogram, default 300
    """
    xs = list(range(-int(spanmax/2),int(spanmax/2)))
    peakwave = {i:0 for i in xs}
    f = pyBigWig.open(bw,'r'); genome = f.chroms(); subsets = []; peaksall = []
    for chrom in genome: #process each chrom separately. 
        wigs = f.intervals(chrom); peaks = [w[0] for w in wigs]
        subset = [peaks[i:i+nsubset] for i in range(len(peaks)-nsubset)]
        if tight=="tight": subsets = [s for s in subsets if s[-2]-s[1]<=spanmax]
        subsets.extend(subset); peaksall.extend(peaks)
    for s in subsets:
        for i in s:
            k = i-s[int(nsubset/2)]
            if k in peakwave: peakwave[k] += 1
    peakwave[0] = 0; ys = [peakwave[i] for i in xs]; plt.plot(xs,ys,ls='-',lw=1)
    plt.savefig(bw[:-3]+"_"+tight+'_peakperiod.pdf'); plt.close()
    f2 = open(bw[:-3]+"_"+tight+'_peakperiod.txt','w')
    for i in range(len(xs)): f2.write('\t'.join([str(xs[i]),str(ys[i])])+'\n')
    f2.write("Number of windows: " + str(len(subsets))+'\n'); f2.close()
    return peaksall #across all chromosomes
################################################################################





################################################################################
##### 2. plot histograms for inter-peak distances, using cntrs from the bw file
def interpeakdist(peaks,tight):
    figname = bw[:-3]+"_"+tight+'_interpeaks.pdf'
    interpeaks = [peaks[i+1]-peaks[i] for i in range(len(peaks)-1)]
    fig, ax = plt.subplots()
    n,bins,patches = plt.hist(interpeaks,bins=200,range=(0,200))
    plt.savefig(figname)
    f = open(bw[:-3]+"_"+tight+'_interpeaks.txt','w')
    for i in range(len(n)): f.write('\t'.join([str(n[i]),str(bins[i])])+'\n')
    f.close()
    return
################################################################################





################################################################################
##### process bw file, export two pdfs and two txt files of distributions. 
if __name__ == '__main__':
    import sys, pyBigWig
    from matplotlib import pyplot as plt

    if len(sys.argv) < 3: sys.exit("python peakperiodicity.py cntrs.bw tight")
    bw = sys.argv[1]; tight = sys.argv[2] #"tight" or "all". 
    #tight means middle 9 peaks or 8 intervals measure <=400nts.
    #all means all peaks even if they are far apart. 

    # step 1. plot and export periodicity in 400nt windows
    peaks = peakperiod(bw,11,400,tight)
    # step 2. plot interpeak distance distribution
    interpeakdist(peaks,tight)
################################################################################





