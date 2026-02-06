"""
phaseregularity.py, zhipengluchina@gmail.com, 2025-07-11. See also phasing.py
calculates loose peaks with prominences, and regularity score, from positioning.
set parameters to find >=5 peaks in 180nt windows, ~ one cycle of tetramer wrap
find all peaks and rank by prominences, and take the top 5 in each window.

~35min for each SHARCLIP sample. 
essential for subsequent analysis of icSHAPE and eCLIP phasing.
as well as relative location to functional motifs, e.g. splicing 

##### Earlier steps before this phasing regularity calculation:
1. extract gap1 alignments and filter out artifacts from all alignments. 
2. convert gap1 reads to short/long/all mid points
3. calculate density and positioning scores, from which regularity is computed.
##### TO DO LIST:
1. plot the score distribution for different regions, genes, datasets, etc.
2. introns vs. exons? different datasets? Other RBPs may have lower regularity.
3. RBPs that bind mostly mature RNAs should have low regularity. 
4. rRNAs should have low regularity.

##### math of the regularity score. 
periodicity = 1-np.std(ds)/np.mean(ds) #ds, inter-peak distances
geoprom = stats.gmean(ps) #ps, peak prominences
regularity = periodicity**4 * geoprom, in range of [0,1]
the 4th power spreads the periodicity score to the full range of [0,1]
the geometric mean is biased towards smaller values, further spreading scores.

##### exports four files:
1. *loosepeaks.bedgraph
2. *loosepeaks_regularity.bw
3. *loosepeaks_regularity_hist.pdf #each inter-peak interval as 1 data point
4. *loosepeaks_regularity_hist.txt

##### Example on HNRNPC SHARCLIP dataset:
script=~/Dropbox/_scripts/phaseregularity.py;
bg2bw=bedGraphToBigWig
genome=~/Documents/lulab/projects/hg38/hg38_SHARCLIP_genome.txt
pre=HepG2_HNRNPC_asso_gap1_sorted_short
pre=HepG2_HNRNPC_indep_gap1_sorted_short
pre=HepG2_HNRNPC_asso_gap1_PICALM_short
python $script ${pre}positioning.bw ${pre}_loosepeaks.bedgraph $bg2bw $genome &

##### Example on server, using xphaseregularity.slurm, with these variables. 
fo=/project/jianhuib_673/sharclip/sharclip_merged_gapped_reads/gap1/positioning/
file=K562_HNRNPC_asso_gap1_sorted_shortpositioning.bw
bg2bw=/project/zhipengl_72/labshare/kentutils/bedGraphToBigWig
genome=/project/zhipengl_72/labshare/hg38/hg38_SHARCLIP_sizes.txt

##### Example on HNRNPC eCLIP datasets: 
Regularity is much lower for HNRNPC eCLIP, confirming high SHARCLIP efficiency.
cd /Users/lu/Desktop/eCLIP_desktop
script=~/Dropbox/_scripts/phaseregularity.py; pre=HepG2_HNRNPC
genome=~/Documents/lulab/projects/hg38/hg38.chrom456.sizes
python $script ${pre}_1_ACTB_positioning.bw \
${pre}_1_ACTB_loosepeaks.bedgraph $genome &
python $script ${pre}_0_positioning.bw ${pre}_0_loosepeaks.bedgraph $genome &
python $script ${pre}_1_positioning.bw ${pre}_1_loosepeaks.bedgraph $genome &
python $script ${pre}_2_positioning.bw ${pre}_2_loosepeaks.bedgraph $genome &


"""



################################################################################
##### 1. calculate regularity from a set of peaks in a window. 
def computeregularity(peaks,npeaks):
    """
    formula: regularity = geoprom * periodicity**4
    regularity=0 if peak number is smaller than the predefined value. 
    input peaks: [[chrom,pos,prom],...]
    also returns all positions for which regularity is calculated. 
    """
    peaks1 = sorted(peaks,key=lambda x:x[2],reverse=True) #sort by prominence
    peaks2 = sorted(peaks1[:min(len(peaks1),npeaks)],key=lambda x:x[1]) #by pos
    ps = [i[2] for i in peaks2] #prominences
    ds = [peaks2[i+1][1]-peaks2[i][1] for i in range(len(peaks2)-1)] #distances
    geoprom = stats.gmean(ps) #geometric mean
    periodicity = 1-np.std(ds)/np.mean(ds) #1 - coefficient of variation
    regularity = geoprom * periodicity**4
    positions = [i[1] for i in peaks if peaks2[0][1]<=i[1]<= peaks2[-1][1]]
    if len(peaks2) < npeaks: regularity = 0
    return positions, regularity 
################################################################################




    
################################################################################
##### 2. process an input peaks bedgraph file from the findpeaks function. 
def arrayregularity(peaksbg,regbg,winsize,npeaks):
    """
    #multiple regularity of inter-peak distances with mean diff80 scores
    #returns a phasing file.
    winsize = 180, window to evaluate the regularity
    npeaks = 5, number of top peaks to pick for evaluation
    """    
    f = open(peaksbg,'r'); peaks = [] #[[chrom,pos,prominence],...]
    chrom0 = '' #process one chrom at a time
    outf = open(regbg,'w')
    for line in f:
        chrom,pos,_,prom = line.split(); pos,prom = int(pos),float(prom)
        if chrom == chrom0:
            if pos - peaks[-1][1] < winsize: 
                peaks.append([chrom,pos,prom])
                while peaks[-1][1]-peaks[0][1] >= winsize: 
                    positions,regularity = computeregularity(peaks,npeaks)
                    for i in range(len(positions)-1):
                        out = [chrom,positions[i],positions[i+1],regularity]
                        outf.write('\t'.join([str(j) for j in out])+'\n')
                    peaks = peaks[1:] #remove first element
            else: peaks = [[chrom,pos,prom]] #next pos too far.
        else: chrom0 = chrom; peaks = [[chrom,pos,prom]]
    f.close(); outf.close()
    return
################################################################################





################################################################################
##### 3. merge lines with identical coordinates, e.g. taking max of values. 
def mergeid(infile,outfile,operation): 
    inf,outf = open(infile, 'r'), open(outfile, 'w')
    prevl,prevc,prevv = '',[],[] #previous line, coord and value
    scores = [] #list of scores for plotting histograms
    for line in inf:
        record = line.split()
        if int(record[1]) > int(record[2]): continue #temp solution for a bug
        if record[:3] == prevc: prevv.append(float(record[3]))
        else:
            if not prevl:
                prevl,prevc,prevv = line,record[:3],[float(record[3])]
                continue
            if operation == "min": v = min(prevv)
            elif operation == "max": v = max(prevv)
            elif operation == "sum": v = sum(prevv)
            elif operation == "mean": v = sum(prevv)/len(prevv)
            elif operation == "median": v = np.median(prevv)
            else: sys.exit("Allowed operations: min,max,sum,mean,median.")
            outf.write('\t'.join(prevc+[str(v)])+'\n'); scores.append(v)
            prevl,prevc,prevv = line,record[:3],[float(record[3])]

    #process last group of lines
    if operation == "min": v = min(prevv)
    elif operation == "max": v = max(prevv)
    elif operation == "sum": v = sum(prevv)
    elif operation == "mean": v = sum(prevv)/len(prevv)
    elif operation == "median": v = np.median(prevv)
    else: sys.exit("Allowed operations: min,max,sum,mean,median.")
    outf.write('\t'.join(prevc+[str(v)])+'\n'); scores.append(v)
    inf.close(); outf.close()
    return scores
################################################################################





################################################################################
##### 4. plot regularity scores in a histogram. Different eCLIP vs. SHARCLIP
def plotreghist(scores,regbg):
    scores = [i for i in scores if i] #remove the zeros
    nbins = 100
    edges = [i/nbins for i in range(nbins+1)]
    figname = regbg[:-9] + '_hist.pdf'
    n,bins,patches = plt.hist(scores,bins=edges,density=True); plt.close()
    n = [i/sum(n) for i in n]
    plt.plot(bins,n+[0])
    plt.savefig(figname); plt.close()
    listfile = regbg[:-9] + '_hist.txt'
    f = open(listfile,'w')
    for i in range(len(n)): f.write('\t'.join([str(bins[i]),str(n[i])])+'\n')
    f.close()
    return 
################################################################################
 




################################################################################
##### 5. merge consecutive records with identical values
def mergeconsec(infile,outfile):
    #input and output: bedgraph files
    f1 = open(infile,'r'); f2 = open(outfile,'w')
    record = f1.readline().split(); coord0,v0 = record[:3],record[3] #previous
    for line in f1:
        record = line.split(); coord,v = record[:3],record[3]
        if coord[1] == coord0[2] and v0 == v: coord0[2] = coord[2] #update
        else: f2.write('\t'.join(coord0+[v0])+'\n'); coord0,v0 = coord,v
    f2.write('\t'.join(coord0+[v0])+'\n')
    f1.close(); f2.close()
    return 
################################################################################





#ignore. Replaced by the regularity method. 
################################################################################
##### 5. oligomer/tetramer array, instead of monomer, positioining.
def getdiff80(pscores,monosize,chrom,start,end):
    """
    compute diff of top/bottom 10 percentile values over 50nt sliding windows.
    pscores: positioning score in a region. 
    output 1: bed file of diff 90% - 10% values, in monosize windows, centered
    output 2: annotation of diff strength as a bedgraph file
    weak: [0,0.2] intron light gray
    medium: [0.2,0.4], thin dark gray
    strong: [0.4,1], thick black
    """
    monosize = 50
    outf = open("PICALM_diff80.bedgraph", 'w')
    outlist = [[chrom,0,0,0]] #[chrom,start,end,value] for bedgraph
    for i in range(len(pscores)-monosize):
        v = sorted(pscores[i:i+monosize])
        if v[0]<0: continue
        c1 = int(i+monosize/2)+start #center of the window
        diff80 = [chrom,c1,c1+1,v[int(monosize*0.9)]-v[int(monosize*0.1)]]
        if diff80[1] == outlist[-1][2]:
            if diff80[3] == outlist[-1][3]: outlist[-1][2] += 1
            else: outlist.append(diff80)
        else: outlist.append(diff80)
    for i in outlist[1:]: outf.write('\t'.join([str(j) for j in i])+'\n')
    outf.close()

    bedf = open("PICALM_diff80_anno.bedgraph", 'w')
    #chrom,start,end,name,score,strand,tstart,tend,rgb,bcount,bstarts,bsizes
    #use thickstart,thickend. thick black, thin dark gray, intron light gray
    bedlist = [[chrom,0,0,0]] #[[chrom,start,end,note]]
    for i in outlist[1:]: #convert to 4 types: blank, weak, medium, strong
        x=0 if i[3]==0 else(0.2 if 0<i[3]<0.2 else(0.4 if 0.2<=i[3]<0.4 else 1))
        if i[1]==bedlist[-1][2] and x==bedlist[-1][3]: bedlist[-1][2] = i[2]
        else: bedlist.append(i[:3]+[x])
    for i in bedlist: bedf.write('\t'.join([str(j) for j in i])+'\n')
    bedf.close()    
    return 
################################################################################






################################################################################
##### 6. process input positioning.bw file
if __name__ == "__main__":

    import os; os.environ['OPENBLAS_NUM_THREADS'] = '1' #override a bug on CARC
    import sys, datetime
    import numpy as np
    import subprocess as sp
    from matplotlib import pyplot as plt
    from scipy import stats
    from struclib import peakfinder


    if len(sys.argv) < 4:
        sys.exit("python phaseregularity.py positionbw peaksbg bg2bw genome")
    positionbw,peaksbg,bg2bw,gfile = sys.argv[1:5]
    chunk = 1000; winsize = 180; npeaks = 5;
    uniq = 5 #>5 values in each 20nt window

    #A. identify loose peaks from positioning bigwig. Relaxed standards.
    print(str(datetime.datetime.today())+"\tFinding loosely defined peaks")
    params={"height":0.2,"distance":30,"prominence":0.05,
            "width":(5,40),"wlen":100} #standard parameters for finding peaks
    params={"height":0.1,"distance":10,"prominence":0.05,
            "width":(5,40),"wlen":60} #finding loose peaks to compute regularity
    peakfinder(positionbw,chunk,peaksbg,params,uniq)

    #B. sort the peaks bedgraph file and merge identical peaks
    print(str(datetime.datetime.today())+"\tSorting peaks bedgraph file")
    cmd = "LC_COLLATE=C sort -k1,1 -k2,2n {} > {}".format(peaksbg,peaksbg+'1')
    sp.run(cmd, shell=True)
    sp.run("mv {} {}".format(peaksbg+'1', peaksbg), shell=True)
    mergeid(peaksbg,peaksbg+'1','max')
    sp.run("mv {} {}".format(peaksbg+'1', peaksbg), shell=True)
    
    #C. compute regularity for entire file. Key step in this script. 
    #export *short_loosepeaks_regularity.bedgraph
    print(str(datetime.datetime.today())+"\tCalculating regularity")
    regbg = peaksbg[:-9] + '_regularity.bedgraph'
    arrayregularity(peaksbg,regbg,winsize,npeaks)

    #D. merge identical intervals in the regularity bedgraph, taking max value
    #rename file name to orginal after merging intervals
    print(str(datetime.datetime.today())+"\tMerging identical intervals")
    cmd = "LC_COLLATE=C sort -k1,1 -k2,2n {} > {}".format(regbg,regbg+'1')
    sp.run(cmd, shell=True)
    sp.run("mv {} {}".format(regbg+'1', regbg), shell=True)
    scores = mergeid(regbg,regbg+'1','max')
    sp.run("mv {} {}".format(regbg+'1',regbg), shell=True)
    
    #E. plot distribution of regularity scores in a histogram/line plot
    #input scores from last step, output histogram pdf file and txt file.
    print(str(datetime.datetime.today())+"\tPlotting regularity scores")
    plotreghist(scores,regbg)

    #F. merge consecutive intervals with identical values.
    #rename file name to orginal. Only do this after plotting the histogram
    print(str(datetime.datetime.today())+"\tMerging consecutive intervals")
    mergeconsec(regbg,regbg+'1')
    sp.run("mv {} {}".format(regbg+'1', regbg), shell=True)
    
    #G. convert bedgraph to bigwig. Use external command
    print(str(datetime.datetime.today())+"\tConverting score bedgraph to bw")
    cmd = "LC_COLLATE=C sort -k1,1 -k2,2n {} > {}".format(regbg,regbg+'1')
    sp.run(cmd, shell=True)
    sp.run("mv {} {}".format(regbg+'1', regbg), shell=True)
    cmd = "{} {} {} {}".format(bg2bw,regbg,gfile,regbg[:-7]+'w')
    sp.run(cmd, shell=True); sp.run("rm {}".format(regbg), shell=True)

    #H. convert loose peaks bedgraph to bigwig
    print(str(datetime.datetime.today())+"\tConverting peaks bedgraph to bw")
    cmd = "{} {} {} {}".format(bg2bw,peaksbg,gfile,peaksbg[:-7]+'w')
    sp.run(cmd, shell=True); sp.run("rm {}".format(peaksbg), shell=True)
################################################################################















