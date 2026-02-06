"""
spliceswitchscan.py, zhipengluchina@gmail.com, 2025-10-22.
see spliceswitch.py and spliceblockSE.py for related functions.
Abundance and structure stability are important ranking metrics.

##### Basic design: use non-overlapping 30nt windows to identify blockers/loops.
On the exon side, we pick SS5+exon+SS3X to increase statistical power. 
We will start with two-MXE groups and then extend to others later.
Multi-MXE clusters require special treatment, but most of them are not real.  
Instead of exporting or plotting DGs, here we export and plot windows.
raw reads are plotted without normalization. Conservations are plotted side by
side, naturally. We can also separate the SS5 and SS3X later.

How do we find alt. loops: use dominant DGs (e.g. top3) in the central intron
and project to both sides? Usually we will be able to find very few real ones.
Together, these two approaches will identify most of the common mechanisms. The
alt loops analysis can use the TCF7L2 example as a model. 

Another note: simple blockers for SE can be presented in a similar way. The
exported bedgraph files can be easily rearranged for pretty display. 

##### As a reference: basic procedure from spliceblockSE.py. 
#A. get gap1 cov in the two introns as norm standards.
#B. get all DGs overlapping SS5, a single nt position
#C. get models for top DGs that overlap SS5
#D. get conservation for top duplexes (not DGs) overlapping SS5
#E. get all DGs overlapping BP+PPT+SS3, e.g. a 40nt region
#F. get models for top DGs that overlap SS3X
#G. get conservation for top duplexes (not DGs) overlapping SS3X

##### Steps, following the spliceblockSE.py approach:
1. no need to get intron coverage for normalization
2. windows blocking the group of MXEs, including SS5 and SS3X.
3. build structure models,
4. also need the conservation as additional filter. focus on the intronic. 
5. load files into IGV for visualization. 
6. Results can be plotted in two bedgraphs, one for SS3X, one for SS5. 

Output different from the standard format, given this design.
MFE can be exported as a bed file and displayed side by side. 
Each cell type requires a separate bedgraph file for the data.
If we do not use bigbed, we can plot positive and negative values together,
For the 14 HNRNPC asso/indep datasets, we will produce 28 files for blocker cov,
SS5 and SS3X separately, and 4 files for MFE, SS3X and SS5 of two MXEs.

##### Example command: 
cd ~/Documents/lulab/projects/sharclip/; mamba activate crssant
script=~/Dropbox/_scripts/spliceswitchscan.py
MXEanno=~/Documents/grants/2026/_2025_SHARCLIP/Hatje2017/Hatje2017_s2_hg38.bed 
fasta=~/Documents/lulab/projects/hg38/fasta/chr19.fa
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_sorted
DGdiv=${pre}_Birch_T20_sorted_DGdiv.bb
phastCons=~/Documents/lulab/projects/hg38/conservation/hg38.phastCons100way.bw
time python $script $MXEanno $fasta $DGdiv $phastCons

The TPM1 first cluster is not in the Hatje2017 annotation.


##### simple quantification of blocker strengths in a range:
#e.g. chr19 38728801 38728831 -4 ACTN4
cd /Users/lu/Documents/lulab/projects/sharclip

#get all the structure coverage in this region
for f in switchscan*; \
do (awk '$4>0 {sum += $4*($3-$2)} END {print sum}' $f); done
#alternatively, get only structure reads in the anchor region:
for f in switchscan*; \
do (awk '$2>=38710300 && $3<=38710500 && $4<0 {sum += $4*($3-$2)} END {print sum}' $f); done


##### To better calculate ACTN4 E8 inclusion levels, I used these custom
commands directly on the bigwig files:
cd /Users/lu/Documents/lulab/projects/sharclip/_iPSsplice/_bws
for f in *bw; do (bigWigAverageOverBed $f ACTN4_E8.bed ACTN4_$f.txt ); done
cd /Users/lu/Documents/lulab/projects/CLIP/encode
for f in *all.bw; do (bigWigAverageOverBed $f ACTN4_E8.bed ACTN4_$f.txt ); done

"""


################################################################################
def switchscan(clusters1,fadict,DGdiv,phastCons,
               samples,indices,flank,loopmax,win):
    """
    get DGs, make bp models, extract conservation. fadict/phastCons not used yet
    clusters1 = {(name,clusterID):[[SE,asmdict[SE],score],...]}
    asmdict[AScoord]['ASdict'] = {cell:[incm,incs,sum_IJC_SJC],...}
    Pairwise comparison in bedgraphs is only compatible with two-MXE clusters. 
    It is important to use whole exon plus SS3X,SS5
    so that the coverage is enough, based on the example TCF3.
    flank=40, extending from SS3 to include PPT and BP, i.e. SS3X region.
    loopmax=50, not used when finding blockers on entire exons. 
    win=30, non-overlapping window size, takes ~1000 steps for 30kb region.
    """
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'; blockdict,loopdict = {},{}
    DGdivf = pyBigWig.open(DGdiv,'r'); phastf = pyBigWig.open(phastCons,'r')
    def divcomma(string): return [int(i) for i in string.split(',')]
    for k in clusters1:
        
        #A. identify minimal introns and SS3X/SS5 in the MXE cluster. 
        name,cluster = k; SS3Xs,SS5s = [],[]
        SEs = [i[0] for i in clusters1[k]]; chrom,strand = SEs[0][:2]
        exons = [[list(SE[2:4]),list(SE[4:6]),list(SE[6:8])] for SE in SEs]
        exons = merge([j for e in exons for j in e])
        introns = [[exons[i][1],exons[i+1][0]] for i in range(len(exons)-1)]
        MXEs,MXEe = exons[0][0],exons[-1][1] #start/end of cluster
        I = introns; w = win; c = chrom #temp variables for windows. 
        for e in exons[1:-1]:
            ES,EE = e[0],e[1]; SS5 = (c,EE-1,EE) if strand=='+' else (c,ES,ES+1)
            SS3X = (c,ES-flank,ES) if strand=='+' else (c,EE,EE+flank)
            SS5s.append(SS5); SS3Xs.append(SS3X)
        windows = [[(c,j,min(i[1],j+w)) for j in range(i[0],i[1],w)] for i in I]
        #print(MXEs,MXEe,SS5s,SS3Xs)        
        
        #X. testing on single examples: 
        #if chrom!='chr19'or MXEs>1611600 or MXEe<1611500: continue #TCF3
        if chrom!='chr19'or MXEs>38711500 or MXEe<38709500: continue #ACTN4
        #if chrom!='chr19'or MXEs>38728900 or MXEe<38727200: continue #ACTN4
        #if chrom!='chr19'or MXEs>1652615 or MXEe<1609292: continue #PKM
        #if chrom!='chr10'or MXEs>92479200 or MXEe<92478000: continue #IDE
        #if chrom!='chr10'or MXEs>121900400 or MXEe<121898400: continue #ATE1
        #if chrom!='chr8'or MXEs>22415000 or MXEe<22408000: continue #SLC39A14
        #if chrom!='chr15'or MXEs>63062300 or MXEe<63060800: continue #TPM1 2
        #if chrom!='chr15'or MXEs>63057500 or MXEe<63042500: continue #TPM1 1
        #if chrom!='chr9'or MXEs>35685400 or MXEe<35684400: continue #TPM2
        #if chrom!='chrX'or MXEs>130157000 or MXEe<130154000: continue #AIFM1
        #if chrom!='chr12'or MXEs>98598000 or MXEe<98594000: continue #AIFM1
        #if chrom!='chr11'or MXEs>73721000 or MXEe<73716000: continue #RAB6A
        
        #if chrom!='chr12'or MXEs>2513000 or MXEe<2493000: continue 


        #B. iterate through defined introns, get DGs blocking/looping MXEs.
        #blockdict[(name,clusterID)][window][exon_index][sample] = coverage
        #loopdict[(name,clusterID)][window][exon_index][sample] = coverage
        #loopdict hard to determine the two regions. See TCF7L2 example 
        r = range(len(SEs)); windows = {j for i in windows for j in i} #flatten
        blockdict[k] = {j:{m:{s:0 for s in samples}for m in r}for j in windows}
        loopdict[k] = {j:{m:{s:0 for s in samples}for m in r}for j in windows}
        
        #C. finding blockers
        entries = DGdivf.entries(chrom,MXEs,MXEe)
        if not entries: continue
        for window in blockdict[k]:
            for entry in entries: #DGs
                s1,e1,string = entry; data = string.split()
                DGID,sizes,starts = data[0],divcomma(data[7]),divcomma(data[8])
                #loopmax only for SS5. SS3X is long, no need for it.
                if overlap(s1,s1+sizes[0],*window[1:]) or \
                   overlap(e1-sizes[1],e1,*window[1:]): #>=1 arms overlap window
                    for i in range(len(SS3Xs)): #>=1 arms overlap SS3X/SS5
                        t1 = overlap(s1,s1+sizes[0],*SS5s[i][1:])  #SS5
                        t2 = overlap(e1-sizes[1],e1,*SS5s[i][1:])  #SS5
                        t3 = overlap(s1,s1+sizes[0],*SS3Xs[i][1:]) #SS3X
                        t4 = overlap(e1-sizes[1],e1,*SS3Xs[i][1:]) #SS3X
                        t5 = overlap(s1,s1+sizes[0],*exons[i+1]) #exon
                        t6 = overlap(e1-sizes[1],e1,*exons[i+1]) #exon                        t5 = 
                        t7 = s1+sizes[0]<SS5s[i][1]<e1-sizes[1] and \
                             e1-sizes[1]-(s1+sizes[0])<=loopmax 
                        #if any([t1,t2,t3,t4,t7]): #block SS3X/SS5 and loop
                        if any([t1,t2,t3,t4,t5,t6]): #block SS3X/SS5 and exon
                            for s in indices:
                                v = int(data[indices[s]])
                                blockdict[k][window][i][s] += v
        
        #E. build bp models for top DGs overlapping SS3X/SS5. No need for cons.
        #use longer regions to build models. up to 100nts or longer. 
        for k in blockdict: pass #try this later if necessary.
    print(blockdict)
    return blockdict
################################################################################





################################################################################
"""
should work well for typical clusters of 2 MXEs, but need special tweaking for
more complex structures. 
"""
def anchorscan(DGdiv,samples,indices,win): 
    #A. manually define anchor regions
    #Heatmaps will work but not as clear. Need to plot multiple pairs.                

    #TCF7L2:
    chrom,strand = 'chr10','-'
    anchor,MXEs,MXEe = [113160050,113160250],113158000,113168000

    #DLG3:
    chrom,strand = 'chrX','+'
    anchor,MXEs,MXEe = [70495500,70497150],70492596,70498520





    #B. finding loops, using defined anchors
    windows = [(chrom,i,i+win) for i in range(MXEs,MXEe,win)]
    loopdict = {j:{s:0 for s in samples} for j in windows} #simpler than above
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'; exps = asso,indep
    DGdivf = pyBigWig.open(DGdiv,'r')
    def divcomma(string): return [int(i) for i in string.split(',')]
    entries = DGdivf.entries(chrom,*anchor)
    for window in loopdict:
        for entry in entries:
            s1,e1,string = entry; data = string.split()
            DGID,sizes,starts = data[0],divcomma(data[7]),divcomma(data[8])
            if overlap(s1,s1+sizes[0],*anchor) and \
               overlap(e1-sizes[1],e1,*window[1:]) or \
               overlap(e1-sizes[1],e1,*anchor) and \
               overlap(s1,s1+sizes[0],*window[1:]):
                for s in indices: loopdict[window][s] += int(data[indices[s]])

    #C. export the loops
    fdict = {}
    for c in all7lines:
        for e in exps:
            fdict[c+e] = open("switchscan_anchor_{}.bedgraph".format(c+e),'w')
    for window in loopdict:
        for s in loopdict[window]:
            if s in fdict:
                l = window+(loopdict[window][s],)
                fdict[s].write('\t'.join([str(j) for j in l])+'\n')
    for k in fdict: fdict[k].close()
    return
################################################################################


  


################################################################################
##### works for both alt blockers and loops.
def exportswitches(blockdict,all7lines):
    #blockdict[(name,clusterID)][window][exon_index][sample] = cov
    #clusterID: int; window: (chrom,start,end); exon_index: 0,1,...;
    #sample: e.g. astroHNRNPCindepmul, focus on the 14 HNRNPC experiments
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'; fdict = {} #outfiles
    for c in all7lines:
        for exp in [asso,indep]:
            k=c+exp; fdict[k] = open("switchscan_{}.bedgraph".format(k),'w')
    for k in blockdict:
        for window in blockdict[k]: #window: chrom,start,end
            for i in blockdict[k][window]:
                if i>1: continue #exon_index, focus on first two for now. 
                sign = 1 if i==0 else -1; name,clusterID = k
                for s in blockdict[k][window][i]:
                    if s in fdict:
                        l = window+(blockdict[k][window][i][s]*sign,name)
                        fdict[s].write('\t'.join([str(j) for j in l])+'\n')
    for k in fdict: fdict[k].close()
    return 
################################################################################





################################################################################
##### process input files, merge several functions
if __name__ == '__main__':
    import sys, pyBigWig; import numpy as np
    from datetime import datetime
    from matplotlib import pyplot as plt
    from struclib import readrMATS,readfa,ssmodel,db2pairs,overlap,merge
    from struclib import makesamples,cwd,rmatsfolders,samples,indices,all7lines
    from spliceswitch import filterMXE
    
    if len(sys.argv) < 5:
        params = ()
        sys.exit("python spliceswitchscan.py MXEanno fasta DGdiv phastCons")
    MXEanno,fasta,DGdiv,phastCons = sys.argv[1:5]
    mininc,maxinc = 0,1
    win = 30 #non-overlapping sliding window size to analyze blockers/loops
    flank = 40 #50? extending the 3' splice sites to include BP and PPT.
    loopmax = 50 #max len of loops in DGs and duplexes that act as blockers
    covmin = 0.02 #relative to max of the region? 
    rmatsfolders = ["astro_vs_neuron_outBAM_igvgtf_rmats"]
    #rmatsfolders = ["astro_vs_neuron_example"]


    #3. scan for loops connected to one anchor region:
    #it works for TCF7L2, even though variation is not super strong.
    #find another better example. 
    print(str(datetime.today())+"starting anchorscan")
    #anchorscan(DGdiv,samples,indices,win); sys.exit()


    #1. read fasta and ... data
    fadict = readfa(fasta)

    #1. read MXE data from rMATS output, using SE files, output asmdict. 
    #output asmdict: key=(chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn) #coord.
    #values = {sample:[inclevel,std,sum_IJC_SJC]}.
    print(str(datetime.today())+'\tStarting analysis')
    asmdict = readrMATS(cwd,rmatsfolders,'SE',mininc,maxinc);
    print(str(datetime.today())+'\tFinished reading rMATS data')
    print("Total SEs in 7 lines:", len(asmdict)) #7 lines:241387. 
    #for k in asmdict:
    #    if k[0]=='chr19' and k[2]>=38710000 and k[3]<=38729000: print(k)
    
    #2. read MXE annotations from Hatje 2017, and filter MXE data. 
    #clusters1, {(name,clusterID):[[SE,asmdict[SE],score],...]}
    #test plotting of inclusion levels. Higher variations are more useful. 
    clusters1 = filterMXE(asmdict,MXEanno,phastCons); #plot_MXEstats(clusters1)
    print(str(datetime.today())+'\tFinished filtering MXE') #check example: 

    #3. quantify  MXE annotations from Hatje 2017, and filter MXE data.
    params=(clusters1,fadict,DGdiv,phastCons,samples,indices,flank,loopmax,win)
    blockdict = switchscan(*params)
    exportswitches(blockdict,all7lines)
    
################################################################################
    









