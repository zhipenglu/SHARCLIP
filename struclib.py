"""
struclib.py, zhipengluchina@gmail.com
started collecting on 2025-07-01

A set of commonly used functions for RNA structure / interaction analysis.
CIGARdiv(CIGAR), returns Glen, gaps, gaplens, segs, Mlens, Qlens, Rlens
tointerval(line), returns intvls #[[RNAME, STRAND, LEFT, RIGHT, LEN], ...]
readfa(fa): read a fasta file and return a dict
bam2DGs(sam,gene): extract DGs from bam and return a dict
rc(seq): returns the reverse complement
extractgenes(minI), minIs, {genename: [chrom,strand,[starts],[ends]], ...}
readbed6(bed6), returns genesdict: {name:[chrom,start,end,score,strand], ...}
overlap(al, ar, bl, br), returns overlap
ssmodel(chrom, intvl, strand, fastadict), returns models, []
peakfinder(bw,chunk,outbg), exports a bedgraph of peaks. 
... 
"""



import os; os.environ['OPENBLAS_NUM_THREADS'] = '1' #for CARC
import re, sys, datetime
import numpy as np
from itertools import product




################################################################################
##### 1. process each line to get coordinates
#all types of alignments are considered. 
def CIGARdiv(CIGAR):
    #considers all CIGAR operations [MINDSPH=X]
    #divide a cigar string to N-separated segments, used in functions below
    #return these results:
    #gaps, a list of N strings, e.g. ['100N', '200N']
    #segs, a list of segments, each with all possible operations except N
    #gaplens, a list of N lengths, e.g. [100, 200]
    #Glen, genomic length of the entire alignment (consumed ref) [MDN=X]
    #Mlens, segment lengths of matches only [M=X]
    #Rlens, segment lengths of consumed Reference [MD=X]
    #Qlens, segment lengths of consumed Query [MIS=X]
    #example: 
    #Glen,gaps,gaplens,segs,Mlens,Qlens,Rlens = CIGARdiv(CIGAR)
    #or replace unwanted output with '_' 
    gaps=re.findall(r'\d+N', CIGAR)
    gaplens=[int(gap[:-1]) for gap in gaps] #gap lengths 
    segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
    Rlens=[sum([int(i[:-1]) for i in re.findall(r'\d+[MD=X]',s)])for s in segs]
    Mlens=[sum([int(i[:-1]) for i in re.findall(r'\d+[M=X]',s)])for s in segs]
    Qlens=[sum([int(i[:-1]) for i in re.findall(r'\d+[MIS=X]',s)])for s in segs]
    Glen=sum(gaplens+Rlens)
    return Glen, gaps, gaplens, segs, Mlens, Qlens, Rlens
################################################################################





################################################################################
##### 2. convert sam/bam line to intervals
def tointerval(line): #requires CIGARdiv
    #given an alignment line, make a list of intervals, all intervals included
    align = line.split(); RNAME, POS = align[2], int(align[3])
    STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
    _,gaps,gaplens,segs,_,_,Rlens = CIGARdiv(align[5]); intvls = [] 
    for i in range(len(segs)): #n gaps and n+1 matches
        intvls.append([RNAME, STRAND, POS+sum(gaplens[:i])+sum(Rlens[:i]), \
                       POS+sum(gaplens[:i])+sum(Rlens[:i+1]), Rlens[i]])
    return intvls #[[RNAME, STRAND, LEFT, RIGHT, LEN], ...]
################################################################################






################################################################################
##### 1. calculate Gini index
def gini(x):
    """ Compute Gini coefficient of array of values. Input is numpy array. """
    if min(x) < 0 or max(x) == 0: print("Only accepts non-negatives"); return 0
    d = sum([np.sum(np.abs(xi - x[i:])) for i,xi in enumerate(x[:-1], 1)])
    return d / (len(x)**2 * np.mean(x))
################################################################################





################################################################################
##### 3. read the fasta file, e.g., hg38.fa entire file
def readfa(fa):
    #output dictionary: fadict {chrnome: seq}
    f = open(fa, 'r'); chroms,falist,fadict = [],[],{} 
    for line in f:
        if line[0]=='>': falist.append([]); chroms.append(line.split()[0][1:])
        else: falist[-1].append(line.split()[0])
    f.close()
    for i in range(len(chroms)): fadict[chroms[i]] = ''.join(falist[i])
    return fadict
################################################################################





################################################################################
##### extract DGs from sam files and return a dict. copied from clusterDGs.py 
def bam2DGs(sam,geneinfo):
    """
    input indexed bam file; geneinfo: [chrom,start,end]
    output: DGdict for each gene
    key=DGID, value=[info,{sampleID:count}]; info: [l1],[l2],[r1],[r2]
    """
    import pysam
    samf = pysam.AlignmentFile(sam, "rb")
    DGdict,DGdict1 = {},{}    
    for line in samf.fetch(*geneinfo):
        line = line.tostring(); record = line.split()
        QNAME,DG = record[0],record[-2][5:]
        sample = '-'.join(QNAME.split('-')[1:]) #print(sample)
        if DG not in DGdict: DGdict[DG] = [[],[],[],[],{}]
        #DGdict element: [left1s,left2s,right1s,right2s,{sample:count}]
        if sample not in DGdict[DG][4]: DGdict[DG][4][sample]=0
        try: intvl1, intvl2 = tuple(tointerval(line)) 
        except: continue
        if intvl1[2]<geneinfo[1] or intvl2[3]>geneinfo[2]: continue
        DGdict[DG][0].append(intvl1[2]); DGdict[DG][1].append(intvl1[3])
        DGdict[DG][2].append(intvl2[2]); DGdict[DG][3].append(intvl2[3])
        DGdict[DG][4][sample]+=1
        
        #if sample == "HepG2HNRNPCindepmul":
        #    print(DG, [(i-geneinfo[1])/5 for i in intvl1[2:4]+intvl2[2:4]])
        #[RNAME,STRAND,LEFT,RIGHT,LEN]
    for k in DGdict:
        if not DGdict[k][0]: continue
        medians = [float(np.median(i)) for i in DGdict[k][:4]]
        DGdict1[k] = [medians,DGdict[k][4]]
    return DGdict1
################################################################################





################################################################################
def rc(seq):
    #converts lowercase to upper case, converts RNA to DNA
    rule = {"A":"T", "a":"T", "T":"A", "t":"A", "U":"A", "u":"A", "G":"C",
            "g":"C", "C":"G", "c":"G", "Y":"R", "y": "R", "R":"Y", "r":"Y",
            "N":"N", "n":"N", "-":"-"}
    return "".join([rule[base] for base in reversed(seq)])
def rccs(seq):
    #case sensitive conversion, converts RNA to DNA
    rule = {"A":"T", "a":"t", "T":"A", "t":"a", "U":"A", "u":"a", "G":"C",
            "g":"c", "C":"G", "c":"g", "Y":"R", "y": "r", "R":"Y", "r":"y",
            "N":"N", "n":"n", "-":"-"}
    return "".join([rule[base] for base in reversed(seq)])
################################################################################







################################################################################
##### 2. extract all minI annotations
def readminI(minI):
    #extract minimal introns to make a dictionary
    minIf = open(minI, 'r')
    minIs = {} #key=genename, value=[chrom,strand,[starts],[ends]] 
    for line in minIf:
        chrom,start,end,name,score,strand = line.split()
        if name not in minIs: minIs[name] = [chrom,strand,[],[]]
        minIs[name][2].append(int(start)); minIs[name][3].append(int(end))
    minIf.close()
    return minIs
################################################################################





################################################################################
def readbed6(bed6):
    """
    input BED6: chrom chromStart chromEnd name score strand
    return genesdict: {name:[chrom,start,end,score,strand], ...}.
    """
    genesdict = {}; genecount = 0
    f = open(bed6, 'r')
    for line in f:
        record= line.rstrip().split('\t'); genecount+=1
        chrom,name,strand = record[0],record[3],record[5]
        start,end,score = int(record[1]),int(record[2]),int(record[4])
        genesdict[name] = [chrom,start,end,score,strand]
    f.close()
    print("Number of genes:", genecount)
    return genesdict
################################################################################







################################################################################
##### 2. overlap function
def overlap(al, ar, bl, br):
    #two intervals a [al, ar] and b [bl, br]; returns the ratio of overlap
    #copied from DGcombo.py, changing a variable name
    if ar < bl or al > br: return 0
    else:
        alength, blength = ar - al, br - bl
        coords = sorted([al, ar, bl, br])
        overlaplength = float(coords[2] - coords[1])
        overlapmin = min(overlaplength/alength, overlaplength/blength)
        return overlapmin
################################################################################





################################################################################
##### 2. overlap function
def overlapmax(al, ar, bl, br):
    #two intervals a [al, ar] and b [bl, br]; returns the ratio of overlap
    #copied from DGcombo.py, changing a variable name
    if ar < bl or al > br: return 0
    else:
        alength, blength = ar - al, br - bl
        coords = sorted([al, ar, bl, br])
        overlaplength = float(coords[2] - coords[1])
        overlap1 = max(overlaplength/alength, overlaplength/blength)
        return overlap1
################################################################################






################################################################################
##### 2. copied from DGcombo.py, modified to consider neg strand and output MFE 
def ssmodel(chrom, intvl, strand, fastadict): 
    #use RNAcofold to obtain structure model; intvl: [l1,l2,r1,r2]
    #used rccs() instead of rc(). 
    #structure prediction is performed on the '+' strand, or rc of '-' strand.
    #for '-' strand, the '+' strand seq is returned, with matching dot-brackets

    import subprocess as sp
    i = intvl
    seq1,cons1 = fastadict[chrom][i[0]:i[1]],"<"*(i[1]-i[0])
    seq2,cons2 = fastadict[chrom][i[2]:i[3]],">"*(i[3]-i[2])
    seq,cons = seq1+"&"+seq2,cons1+"&"+cons2
    if strand == '-':
        seq,cons = rccs(seq2)+'&'+rccs(seq1),"<"*(i[3]-i[2])+'&'+">"*(i[1]-i[0])
    #if 'a' in seq or 't' in seq or 'c' in seq or 'g' in seq: return '','',0
    p = sp.Popen(["RNAcofold", "-C"], stdout=sp.PIPE, stdin=sp.PIPE)
    IN = (seq + "\n" + cons).encode()
    helix = p.communicate(input=IN)[0].decode().split('\n')
    struc = helix[1][0:len(helix[0])] #structure: e.g., ...(((....)))...
    if strand == '-':
        struc = struc.replace(')','x').replace('(',')').replace('x','(')[::-1]
    seq = seq1+"&"+seq2 #return '+' strand, even if seq is on the '-' strand
    mfe = float(helix[1][len(helix[0])+2:-1])
    return seq, struc, mfe
    #seq,struc,mfe = ssmodel(chrom,intvl,strand,fastadict)
################################################################################





################################################################################
##### 1. find peaks in bw files using the signal.find_peaks().
"""
detailed notes on scipy.signal.find_peaks().
Alternatively use find_peaks_cwt() to smoothen the noisy signal. 
find_peaks(x, height=None, threshold=None, distance=None, prominence=None,
           width=None, wlen=None, rel_height=0.5, plateau_size=None)
returns ndarray peaks and dict properties: 
‘peak_heights’: absolute values, not so useful
‘left_thresholds’, ‘right_thresholds’: not so useful
‘prominences’, ‘right_bases’, ‘left_bases’: useful for positioning analysis
‘widths’, ‘width_heights’, ‘left_ips’, ‘right_ips’: used to calculate centers
‘plateau_sizes’, left_edges’, ‘right_edges’. Not necessary. 

##### Key considerations:
1. 'x', positioning score, is in range [0,1], reducing bias. 
2. 'height' should be >0.18, because 11/61 ~ 0.18. Alsolute height
3. 'threshold' is difficult to set at the moment
3. 'distance' between neighbor peaks set to 30, lower limit. Not for ips
4. 'prominence' is calculated internally
5. 'width' set based on inspection: (5,30), min/max at half height
6. 'wlen' set based on inspection: 100. Should be enough for even wider valleys
7. 'rel_height' set to default 0.5 to calculate 'left_ips' and 'right_ips'
8. 'peaks' output are often skewed. Instead we use (left_ips+right_ips)/2

peaks,properties=find_peaks(x,distance=30,width=(5,30),wlen=150,rel_height=0.5)
ipsl,ipsr=properties["left_ips"],properties["right_ips"]
modes = [int((ipsl[i]+ipsr[i])/2) for i in range(len(ipsl))]

##### Use properties to find tetramers/triads: strong positioning and phasing 
1. variance of centers, in place of 'peaks', to assess phasing
2. mean/product of 'prominences', instead of gap1 pos score (GPS), positioning
3. other RBPs and their specific binding patterns.
"""
def peakfinder(bw,chunk,outbg,params,uniq):
    """
    input: bigwigs of positioning. Read bw chunk by chunk.
    output: bedgraphs of centers, based on half heights of the mountains
    output format: [chrom, start, end, 1]
    chunk: default ~1000nt range, give or take the first/last peak positions
    Only retain ones where in the 20nt window there are >5 unique values
    Peak output from signal.find_peaks() is likely biased due to skewed shape.
    """
    import pyBigWig
    from scipy import signal
    if not params: #default as follows, otherwise provide custom ones
        params = {"height":0.2,"distance":30,"prominence":0.05,
                  "width":(5,40),"wlen":100}
    bwf = pyBigWig.open(bw); bgf = open(outbg,'w'); chroms = bwf.chroms() 
    for chrom in chroms: #chroms: {name:length}
        print(str(datetime.datetime.today())+"\tProcessing", chrom)
        lastc = int(1E9) #initialize last peak center in the window
        for i in range(0,chroms[chrom],chunk):
            #adjust start and end depending on last chunk
            start,end = min(i,lastc),min(i+chunk,chroms[chrom]) #region in bw
            wig = bwf.intervals(chrom,start,end) #wig: [[start,end,value],...]
            if not wig: continue
            start1,end1 = wig[0][0],wig[-1][1] #region to fill with values
            d = [0 for i in range(0,end1-start1)]
            for w in wig: #w: [start,end,value]
                for j in range(w[0]-start1,w[1]-start1): d[j] = w[2]
            peaks,properties = signal.find_peaks(d,**params)
            ipsl,ipsr = properties["left_ips"],properties["right_ips"]
            proms = properties["prominences"]
            cs = [int((ipsl[i]+ipsr[i])/2)+start1 for i in range(len(ipsl))]
            if not cs: continue
            for i in range(len(cs)): # peak centers
                pos1,pos2 = max(0,cs[i]-10),min(cs[i]+10,chroms[chrom])
                vs = [j[2] for j in bwf.intervals(chrom,pos1,pos2)]
                if len(set(vs))>uniq: #peak is 'good' if n of unique values > 5.
                    out = [chrom,str(cs[i]),str(cs[i]+1),str(proms[i])]
                    bgf.write('\t'.join(out)+'\n')
            lastc = start1+cs[-1] #last peak center position
    bwf.close(); bgf.close()
    return
################################################################################




#discarded. 
################################################################################
##### 1. read the genome size file
def getgenomesize(sizefile):
    f = open(sizefile, 'r'); genome = {}
    for line in f: record = line.split(); genome[record[0]] = int(record[1])
    f.close()
    return genome
################################################################################





################################################################################
##### 1. read the chunks bed file. 
def bed3tolist(bed3):
    f = open(bed3,'r'); List = []
    for line in f: r = line.split(); List.append([r[0],int(r[1]),int(r[2])])
    return List
################################################################################





################################################################################
##### 2. build a tree of intervals for repeats
def bed2trees(bed):
    """
    general purpose convert a bed file to a dict of trees, key=chrom
    strand, if present, divides each chrom to two trees, otherwise use '.'
    other bed fields [3:] are stored as additional information
    """
    import intervaltree as itree
    f = open(bed,'r'); trees = {}
    for line in f:
        record = line.split(); start,end = int(record[1]),int(record[2])
        chrom = record[0]; strand = '.' if len(record) < 6 else record[5]
        info = record[3:]; k = (chrom,strand)
        if k not in trees: trees[k] = itree.IntervalTree()
        trees[k].addi(start,end)
    return trees
################################################################################





################################################################################
def density_scatter(x,y,outfile,xlims,ylims,ax=None,sort=True,bins=20):
    print("Input x and y must be numpy.array, not a list!")
    #https://stackoverflow.com/questions/20105364/\
    #how-can-i-make-a-scatter-plot-colored-by-density/53865762#53865762
    #Scatter plot colored by 2d histogram
    #modified for general usage
    from scipy.interpolate import interpn
    from matplotlib import pyplot, cm, colors
    maxv = max([max(x),max(y)])
    if ax is None: fig, ax = pyplot.subplots()
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density = True )
    z = interpn((0.5*(x_e[1:]+x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])),
                data,np.vstack([x,y]).T, method="splinef2d", bounds_error=False)
    z[np.where(np.isnan(z))] = 0.0 #To be sure to plot all data

    # Sort the points by density, so that the densest points are plotted last
    if sort: idx = z.argsort(); x, y, z = x[idx], y[idx], z[idx]
    ax.scatter(x,y, c=z,)
    norm = colors.Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax)
    cbar.ax.set_ylabel('Density'); maxcov = max(x+y)
    #pyplot.gca().set_aspect("equal")
    pyplot.xlim(*xlims); pyplot.ylim(*ylims); pyplot.savefig(outfile)
    return ax
################################################################################






################################################################################
##### 1. process input sam files, allowing selection of a shorter region
def bam2matrix(sam,s1,s2,radius,
               chrom,start,rnasize,binsize,spanmax,matrixnorm,option):
    """
    Each sam file includes data for multiple samples.
    radius: around midpoint to select for higher resolution, default=7
    s1,s2: sample names, e.g. HepG2HNRNPCassomul
    chrom,start,rnasize: define the region to be analyzed.
    binsize: for plotting purposes. defaults=5,20,100,300, etc. 
    matrix computation is set by option: all, median, or trim (medians +/7nts)    
    """
    import pysam
    f = pysam.AlignmentFile(sam, "rb");
    s,b = start,binsize; ncol = int(rnasize/b+1)
    M0,M1 = {},{} #cov of each arm on diagonal or contacts, {sample:np.array}
    Ms = {} #s for short, cov of contacts with span <= certain distance
    for line in f.fetch(chrom,start,start+rnasize):
        line = line.to_string(); record = line.split()
        sample = '-'.join(record[0].split('-')[1:])
        if s1 not in sample and s2 not in sample: continue
        for M in [M0,M1,Ms]: M[sample] = M.get(sample,np.zeros((ncol, ncol)))
        intvl1,intvl2 = tuple(tointerval(line))
        l1,l2,r1,r2 = *intvl1[2:4],*intvl2[2:4]
        if min([l1,l2,r1,r2])<s or max([l1,l2,r1,r2])>s+rnasize: continue
        ml,mr = int((l1+l2)/2),int((r1+r2)/2) #medians
            
        if option == 'median': # option 1. fast but too pixelated. 
            M0[sample][int((ml-s)/b),int((ml-s)/b)]+=1 #coverage
            M0[sample][int((mr-s)/b),int((mr-s)/b)]+=1 #coverage
            M1[sample][int((ml-s)/b),int((mr-s)/b)]+=1 #contacts
            #does not produce the Ms matrix
        elif option in ['all','trim']:
            # options 2-3. slow but contains all or most information
            if option == 'trim': # median +/- radius, e.g. 7nt, more continuous
                l1,l2,r1,r2 = ml-radius,ml+radius,mr-radius,mr+radius
            range1,range2 = list(range(l1,l2)),list(range(r1,r2)) #range lists
            for (x,y) in list(product(range1,range1))+\
                list(product(range2,range2)):
                M0[sample][int((x-s)/b),int((y-s)/b)] +=1 #coverage
            for (x,y) in product(range1,range2):
                M1[sample][int((x-s)/b),int((y-s)/b)] +=1 #contacts
            if mr-ml>=spanmax: continue
            for (x,y) in product(range1,range2):
                Ms[sample][int((x-s)/b),int((y-s)/b)] +=1 #contacts
        else: sys.exit("Error: set option to 'all', 'median' or 'trim'.")
    f.close()

    if matrixnorm: 
        M1norm,Msnorm = {},{}
        for sample in M1: #normalize by 90th percentile, i.e. matrixnorm=0.9
            if not M1[sample].sum(): continue #ignore samples with all zeros.
            list1 = sorted([i for i in M1[sample].flatten().tolist() if i])
            r = list1[int(len(list1)*matrixnorm)]
            if list1: M1norm[sample],Msnorm[sample] = M1[sample]/r,Ms[sample]/r
        M1,Ms = M1norm, Msnorm
    return M0, M1, Ms
################################################################################






################################################################################
##### 1. process input sam files, allowing selection of a shorter region
#compared to bam2matrix, this version subsamples the dataset with more reads
def bam2matrixSub(sam,s1,s2,radius,chrom,start,rnasize,binsize,spanmax,option):
    """
    Each sam file includes data for multiple samples.
    radius: around midpoint to select for higher resolution, default=7
    s1,s2: sample names, e.g. HepG2HNRNPCassomul
    chrom,start,rnasize: define the region to be analyzed.
    binsize: for plotting purposes. defaults=5,20,100,300, etc. 
    matrix computation is set by option: all, median, or trim (medians +/7nts)
    """
    import random,pysam
    f = pysam.AlignmentFile(sam, "rb");
    s,b = start,binsize; ncol = int(rnasize/b+1)
    M0,M1 = {},{} #cov of each arm on diagonal or contacts, {sample:np.array}
    Ms = {} #s for short, cov of contacts with span <= certain distance

    s1reads,s2reads = [],[] #sample the same number of reads
    for line in f.fetch(chrom,start,start+rnasize):
        line = line.to_string(); record = line.split()
        sample = '-'.join(record[0].split('-')[1:])
        if s1 not in sample and s2 not in sample: continue
        intvls = tuple(tointerval(line))
        if len(intvls)>2: continue
        intvl1,intvl2 = intvls
        l1,l2,r1,r2 = *intvl1[2:4],*intvl2[2:4]
        if min([l1,l2,r1,r2])<s or max([l1,l2,r1,r2])>s+rnasize: continue
        if s1 in sample: s1reads.append(line)
        elif s2 in sample: s2reads.append(line)
    f.close()
    s1count,s2count = len(s1reads),len(s2reads) #only count completely in range
    sless,smore = s1reads,s2reads
    if s1count>s2count: sless,smore = s2reads,s1reads
    reads = sless+random.sample(smore,min(s1count,s2count))

    for line in reads:
        record = line.split(); sample = '-'.join(record[0].split('-')[1:])        
        for M in [M0,M1,Ms]: M[sample] = M.get(sample,np.zeros((ncol, ncol)))
        intvl1,intvl2 = tuple(tointerval(line))
        l1,l2,r1,r2 = *intvl1[2:4],*intvl2[2:4]
        ml,mr = int((l1+l2)/2),int((r1+r2)/2) #medians

        if option == 'median': # option 1. fast but too pixelated. 
            M0[sample][int((ml-s)/b),int((ml-s)/b)]+=1 #coverage
            M0[sample][int((mr-s)/b),int((mr-s)/b)]+=1 #coverage
            M1[sample][int((ml-s)/b),int((mr-s)/b)]+=1 #contacts
            #does not produce the Ms matrix
        elif option in ['all','trim']:
            # options 2-3. slow but contains all or most information
            if option == 'trim': # median +/- radius, e.g. 7nt, more continuous
                l1,l2,r1,r2 = ml-radius,ml+radius,mr-radius,mr+radius
            range1,range2 = list(range(l1,l2)),list(range(r1,r2)) #range lists
            for (x,y) in list(product(range1,range1))+\
                list(product(range2,range2)):
                M0[sample][int((x-s)/b),int((y-s)/b)] +=1 #coverage
            for (x,y) in product(range1,range2):
                M1[sample][int((x-s)/b),int((y-s)/b)] +=1 #contacts
            if mr-ml>=spanmax: continue
            for (x,y) in product(range1,range2):
                Ms[sample][int((x-s)/b),int((y-s)/b)] +=1 #contacts
        else: sys.exit("Error: set option to 'all', 'median' or 'trim'.")
    return M0, M1, Ms
################################################################################





################################################################################
##### 1. process input sam files, allowing selection of a shorter region
def bam2matrix1(sam,s1,radius,chrom,start,rnasize,binsize,spanmax,option):
    """
    Each sam file includes data for multiple samples.
    radius: around midpoint to select for higher resolution, default=7
    s1: sample name, e.g. HepG2HNRNPCassomul
    chrom,start,rnasize: define the region to be analyzed.
    binsize: for plotting purposes. defaults=5,20,100,300, etc. 
    matrix computation is set by option: all, median, or trim (medians +/7nts)
    """
    import pysam
    f = pysam.AlignmentFile(sam, "rb");
    s,b = start,binsize; ncol = int(rnasize/b+1)
    M1 = {} #cov of the contacts, {sample:np.array}
    for line in f.fetch(chrom,start,start+rnasize):
        line = line.to_string(); record = line.split()
        sample = '-'.join(record[0].split('-')[1:])
        if s1 not in sample: continue
        M1[sample] = M1.get(sample,np.zeros((ncol, ncol)))
        try: intvl1,intvl2 = tuple(tointerval(line))
        except: continue
        l1,l2,r1,r2 = *intvl1[2:4],*intvl2[2:4]
        if min([l1,l2,r1,r2])<s or max([l1,l2,r1,r2])>s+rnasize: continue
        ml,mr = int((l1+l2)/2),int((r1+r2)/2) #medians
            
        if option == 'median': # option 1. fast but too pixelated. 
            M1[sample][int((ml-s)/b),int((mr-s)/b)]+=1 #contacts
        elif option in ['all','trim']:
            # options 2-3. slow but contains all or most information
            if option == 'trim': # median +/- radius, e.g. 7nt, more continuous
                l1,l2,r1,r2 = ml-radius,ml+radius,mr-radius,mr+radius
            range1,range2 = list(range(l1,l2)),list(range(r1,r2)) #range lists
            for (x,y) in product(range1,range2):
                M1[sample][int((x-s)/b),int((y-s)/b)] +=1 #contacts
        else: sys.exit("Error: set option to 'all', 'median' or 'trim'.")
    f.close()
    return M1
################################################################################






################################################################################
##### 1. extract AS data from one file. see similar script readrMATS()
#deprecated? 
def extractAS(infile,fields,pmax,incdmin):
    """
    infile: [AS].MATS.JC.txt; len(fields)=23 for A3SS,A5SS,RI,SE; 25 for MXE
    AS types: ['A3SS','A5SS','MXE','RI','SE']
    Build 2-level dict for each file: ASdict[coord][field] = value
    fields: ['IJC1','SJC1','IJC2','SJC2','p','FDR','inc1','inc2','incd']
    coord: 9/11-tuple: (gene,chrom,strand,...), A3SS,A5SS,RI,SE:9, MXE:11
    comparison folders: [list of folders], contains samples compared in rMATS
    value: depends on field, float, or list of floats, or list of ints
    setting pmax=1 and incdmin=0 will keep all the records in the file
    "NA" is converted to 1 for p and FDR, or 0 for incd.
    """
    f = open(infile,'r'); ASdict = {}
    for line in f:
        record = line.split()
        if record[0] == "ID": continue
        x = 0 if len(record) == 23 else 2 #0 for A3SS,A5SS,RI,SE; 2 for MXE
        coord = tuple(record[2:5]+[int(i) for i in record[5:11+x]])
        jcounts = [[int(j) for j in i.split(',')] for i in record[12+x:16+x]]
        incs = [[0 if j=="NA" else float(j) for j in i.split(',')]
                for i in record[20+x:22+x]]
        p = 1 if record[18+x] == 'NA' else float(record[18+x])
        FDR = 1 if record[19+x] == 'NA' else float(record[19+x])
        incd = 0 if record[22+x] == 'NA' else float(record[22+x])
        values = [*jcounts,p,FDR,*incs,incd]
        if p <= pmax and abs(incd) >= incdmin:
            ASdict[coord] = {}
            for i in range(len(fields)):
                ASdict[coord][fields[i]] = values[i]
    f.close()
    print("No. of significant AS events in {}: {}".format(infile,len(ASdict)))
    return ASdict
################################################################################






################################################################################
##### 2. collect AS events from all files and store in dicts
#complicated but worked well. only used in ASmerge.py, not sure if still useful.
def mergeAS(AStypes,folders,fields,pmax,incdmin):
    """
    AStypes: ['A3SS','A5SS','MXE','RI','SE']
    ASdict from extractAS(), ASdict[coord][field] = value
    coord: 9/11-tuple: (gene,chrom,strand,...), A3SS,A5SS,RI,SE:9, MXE:11
    fields: ['IJC1','SJC1','IJC2','SJC2','p','FDR','inc1','inc2','incd']
    """
    #A. Build a 4-level dict: dicts1[AS][folder][coord][field] = value
    cwd = os.getcwd(); dicts1,dicts = {},{}
    for AS in AStypes:
        dicts1[AS] = {}
        for folder in folders: #extract AS events from each folder.
            infile = '{}/{}/{}.MATS.JC.txt'.format(cwd,folder,AS)
            dicts1[AS][folder] = extractAS(infile,fields,pmax,incdmin)

    #B. reorganize to dicts[AS][coord][field][folder] = value
    for AS in AStypes:
        dicts[AS] = {}
        coords = [list(dicts1[AS][folder].keys()) for folder in folders]
        coords = list(set([j for i in coords for j in i])) #flatten list
        for coord in coords:
            dicts[AS][coord] = {}
            for field in fields:
                dicts[AS][coord][field] = {}
                for folder in folders:
                    if coord in dicts1[AS][folder] and \
                       field in dicts1[AS][folder][coord]: 
                        value = dicts1[AS][folder][coord][field]
                        dicts[AS][coord][field][folder] = value
    return dicts #dicts[AS][coord][field][folder] = value
################################################################################







################################################################################
##### this function is different from extractAS(), which processes 1 input file
#will be used to analyze all AS events, replacing extractAS, getSE, etc.

def readrMATS(cwd,folders,AS,mininc,maxinc):
    """
    read splicing data from all rMATS output files
    input: cwd, current working directory,
    input: folders, and a list of rMATS output folders
    infile: [AS].MATS.JC.txt; len(fields)=23 for A3SS,A5SS,RI,SE; 25 for MXE
    AS types: ['A3SS','A5SS','MXE','RI','SE']
    output alternative splicing mechanism dict, or asmdict:
    key = coord, (chrom,strand, and 6/8 numbers) #8/10-tuple
    values = {}, where 'ASdict':{sample:[inclevel,std,sum_IJC_SJC]}.
    later other dicts will be added, primarily DGdict, and related data. 
    standardized variable names will be adopted from now on. 

    A3SS/A5SS 6-11 columns: ESlong EElong ESshort EEshort ESflank EEflank
    RI 6-11 columns: riES riEE ESup EEup ESdn EEdn
    SE 6-11 columns: ES EE ESup EEup ESdn EEdn
    MXE 6-13 columns: ES1 EE1 ES2 EE2 ESup EEup ESdn EEdn

    #for HepG2 and K562, ENCORE coverage was lower in each sample
    #more samples were available, so the results may still be accurate. 
    """
    file = '/{}.MATS.JC.txt'.format(AS); asmdict = {}
    def divcomma(string): return [int(i) for i in string.split(',')]
    for folder in folders:
        f = open(cwd+folder+file,'r'); cell1,_,cell2 = folder.split('_')[:3]
        for line in f:
            record = line.split() #print(cwd+folder+file)
            if 'ID' in record[0] or 'txt' in record[0]: continue 
            x = 0 if len(record) == 23 else 2 #0 for A3SS,A5SS,RI,SE; 2 for MXE
            AScoord = tuple(record[3:5]+[int(i) for i in record[5:11+x]])
            IJC1,SJC1 = divcomma(record[13+x]),divcomma(record[14+x])
            IJC2,SJC2 = divcomma(record[15+x]),divcomma(record[16+x])

            inc1 = [float(i) for i in record[20+x].split(',') if i!='NA']
            inc2 = [float(i) for i in record[21+x].split(',') if i!='NA']
            #if 'NA' in record[20+x] or 'NA' in record[21+x]: continue #no 'NA'
            if not inc1 or not inc2: continue #'NA' allowed in some samples
            inc1m,inc1s = float(np.mean(inc1)),float(np.std(inc1))
            inc2m,inc2s = float(np.mean(inc2)),float(np.std(inc2))
            if AScoord not in asmdict: asmdict[AScoord] = {'ASdict':{}}
            asmdict[AScoord]['ASdict'][cell1] = [inc1m,inc1s,sum(IJC1+SJC1)]
            asmdict[AScoord]['ASdict'][cell2] = [inc2m,inc2s,sum(IJC2+SJC2)]
    asmdict1 = {}
    for AScoord in asmdict:
        ASdict = asmdict[AScoord]['ASdict'];
        inc = [ASdict[s][0] for s in ASdict]
        if mininc<=min(inc)<=maxinc:
            asmdict1[AScoord] = {}
            asmdict1[AScoord]['ASdict'] = asmdict[AScoord]['ASdict']
    return asmdict1 #there will be no filtering if maxinc is set to 1. 
################################################################################





################################################################################
def bedpeto12(bedpe,bed12):
    bedpef = open(bedpe,'r'); bed12f = open(bed12,'w')
    for line in bedpef:
        record = line.strip('\n').split(); chrom = record[0]
        start1, end1 = int(record[1]), int(record[2])
        start2, end2 = int(record[4]), int(record[5])
        if start1==end1 or start2==end2: continue #remove artifacts
        name, score, strand = record[6:9]; RGB = "0,0,0"
        sizes = str(end1-start1)+","+str(end2-start2)
        starts = "0," + str(start2-start1)
        outrecord = [chrom, str(start1), str(end2), name, score, strand, \
                     str(start1), str(start1), RGB, "2", sizes, starts]
        bed12f.write("\t".join(outrecord) + "\n")
    bedpef.close(); bed12f.close(); return 
################################################################################





################################################################################
##### 1. extract bed records from IEIbed file as a dictionary
#chr14 21211306 21211808 HNRNPC  1 - 21211306 21211809 0,0,0 2 99,243 0,260
#194944 records. 
def getIEI(IEIbed): #IEIbed: bed12
    f = open(IEIbed,'r'); IEIdict = {}
    #{(chrom,strand,ES,EE):[bed12,[],[],[],[]]}; ES,EE: skipped exon start/end
    def divcomma(string): return [int(i) for i in string.split(',')]
    for line in f:
        chrom,start,end,name,_,strand,_,_,_,_,sizes,starts = line.split()
        sizes,starts = divcomma(sizes),divcomma(starts)
        ES,EE = int(start)+sizes[0],int(start)+starts[1]
        IEIdict[(chrom,strand,ES,EE)] = [line,[],[],[],[]] #IJC1,SJC1,IJC2,SJC2
    f.close(); return IEIdict 
################################################################################
 




################################################################################
def getEIEIE(IEIbed):
    """
    note: splicecombo.py only makes 4 combinations: minI, maxE, IEI and EIE
    additional ones are made from the above combinations. 
    """
    pass
################################################################################





################################################################################
##### 2. retrieve splicing data from SEfile: IJC1,SJC1,IJC2,SJC2.
#used in spliceloop.py and other related scripts
#only keep junctions with >=minj reads (IJC+SJC)
def getSE(IEIdict,SEfile,minj): #SEfile: 23 fields; minj: min junction counts
    f = open(SEfile,'r')
    def divcomma(string): return [int(i) for i in string.split(',')]
    for line in f: #add splicing data to the IEIdict
        record = line.split()
        if record[0] == "ID": continue
        _,_,name,chrom,strand,ES,EE,_,_,_,_,_,I1,S1,I2,S2,_,_,_,_,_,_,_ = record
        I1sum,S1sum = sum(divcomma(I1)),sum(divcomma(S1))
        I2sum,S2sum = sum(divcomma(I2)),sum(divcomma(S2))
        SE = (chrom,strand,int(ES),int(EE)); sums = [I1sum,S1sum,I2sum,S2sum]
        if SE in IEIdict:
            for i in range(1,5): IEIdict[SE][i].append(sums[i-1])
    f.close(); tmp = {}
    for k in IEIdict:
        n = [sum(i) for i in IEIdict[k][1:]]; sum1,sum2 = sum(n[:2]),sum(n[2:])
        inc1 = 'NA' if sum1<minj else n[0]/sum1
        inc2 = 'NA' if sum2<minj else n[2]/sum2
        if inc1 != 'NA' or inc2 != 'NA': tmp[k] = [IEIdict[k][0],inc1,inc2]
    return tmp #{(chrom,strand,ES,EE):[bed12,inc1,inc2]}
################################################################################






################################################################################
def db2pairs(db):
    #Converts a dot-bracket string to a list of base pair coordinates.
    #Pseudoknots are considered, using these notes: ['()','[]','{}','<>'] 
    #Returns a list of lists, indices are 0-based, e.g. [[],[],[],...]
    #here coordinates are based on the duplex, not the genome
    pairs = []
    for p in ['()','[]','{}','<>']:
        pl,pr = p[0],p[1]; stack = [] # Stores indices of opening parentheses
        for i, char in enumerate(db):
            if char == pl: stack.append(i)
            elif char == pr:
                if stack: openi = stack.pop(); pairs.append([openi, i])
                else: pass #handles exceptions. Not implemented here.
    return pairs
################################################################################




################################################################################
def mergehelices(sublists,i,j):
    #Groups sublists that are chained together head-to-tail.
    #example: sublists = [[(1,9),(2,8)],[(4,6),(3,5)],[(12,18),(13,17)]]
    #Store sublists by their first and last elements for quick lookup.
    #Track visted sublists to prevent reprocessing
    starts,ends = {},{}; result = []; visited = set()
    for sl in sublists: starts[sl[0]] = sl; ends[sl[-1]] = sl
    for sl in sublists:
        if tuple(sl) not in visited:
            current = sl[:]; visited.add(tuple(sl)) #start a new chain
            while (current[-1][0]+i,current[-1][1]-j) in starts: #follow forward
                next_sl = starts[(current[-1][0]+i,current[-1][1]-j)]
                if tuple(next_sl) in visited: break
                current.extend(next_sl); visited.add(tuple(next_sl))
            while (current[0][0]-i,current[0][1]+j) in ends: #follow backward
                prev_sl = ends[(current[0][0]-i,current[0][1]+j)]
                if tuple(prev_sl) in visited: break
                current = prev_sl+current; visited.add(tuple(prev_sl))
            result.append(current)
    return result #same format as the input sublists
################################################################################




################################################################################
##### 1. function to merge overlapping or identical intervals
def merge(intvls):
    #merge overlapping intervals, format: [[x1,y1],...]
    intvls.sort(key=lambda i:i[0]); merged = []
    for i in intvls:
        if not merged or merged[-1][1] < i[0]: merged.append(i)
        else: merged[-1][1] = max(merged[-1][1],i[1])
    return merged
################################################################################






################################################################################
def itermerge(bplist,n): 
    #input bplist: [(x,y),...], does not have to be sorted, uses mergehelices()
    #max internal bulge or loop size, e.g. n = 5
    helices = [[pair] for pair in bplist] #initial helices, one pair per helix
    distances = list(product(range(1,n),range(1,n)))
    distances = sorted(distances,key=lambda x:sum(x))
    for i,j in distances: helices = mergehelices(helices,i,j)
    return helices
################################################################################






        
################################################################################
def altcheck(intvl1,intvl2,tover,talt):
    #returns whether intvls of two DGs overlap, using the overlap() function
    #input intvl1/intvl2, each a 4-element list: [l1,l2,r1,r2]
    #here intvls may be DG arms, or their predicted structures. 
    #tover and talt are defined arbitrarily, and defaults are 0.5 and 0.2
    i11,i12,i13,i14 = intvl1[0],intvl1[1],intvl1[2],intvl1[3]
    i21,i22,i23,i24 = intvl2[0],intvl2[1],intvl2[2],intvl2[3]
    o1,o2 = overlap(i11,i12,i21,i22),overlap(i13,i14,i23,i24)
    o3,o4 = overlap(i11,i12,i23,i24),overlap(i13,i14,i21,i22)
    return (o1>tover)*(o2<talt) or (o1<talt)*(o2>tover) or o3>tover or o4>tover
################################################################################




################################################################################
def altcheckbp(chrom,intvl1,intvl2,fadict,tover,talt):
    #use RNAcofold to obtain structure and test overlap, using overlap()
    #input intvl1/intvl2, a 4-element list: [l1,l2,r1,r2]
    #return a 3-tuple if overlap is above tover
    #maxoverlap: int number of overlapping bps. 
    #ss1/2: dot-bracket format of structure, e.g. ...(((....)))...
    ss,si = [],[] #structures and intervals. si: [[a1,a2,a3,a4],[b1,b2,b3,b4]]
    for intvl in [intvl1,intvl2]: 
        seq,s,mfe = ssmodel(chrom,intvl,strand,fadict); ss.append(s)
        s = s.replace(")", "(").split("&"); sintvl = [0,0,0,0] #why
        #remove dots at begining and end to make the new intervals for brackets
        sintvl[0] = intvl[0] + len(s[0].split("(")[0])
        if s[0][-1] != ".": sintvl[1] = intvl[1]
        else: sintvl[1] = intvl[1] - len(s[0].split("(")[-1])
        if s[1][0] != ".": sintvl[2] = intvl[2]
        else: sintvl[2] = intvl[2] + len(s[1].split("(")[0])
        if s[1][-1] != ".": sintvl[3] = intvl[3]
        else: sintvl[3] = intvl[3] - len(s[1].split("(")[-1])
        si.append(sintvl)
    si11,si12,si13,si14 = si[0][0],si[0][1],si[0][2],si[0][3]
    si21,si22,si23,si24 = si[1][0],si[1][1],si[1][2],si[1][3]
    o1,o2 = overlap(si11,si12,si21,si22),overlap(si13,si14,si23,si24)
    o3,o4 = overlap(si11,si12,si23,si24),overlap(si13,si14,si21,si22)
    maxoverlap = max(o1,o2,o3,o4)
    if maxoverlap > tover: return maxoverlap, ss[0], ss[1]
    else: return 0, ss[0], ss[1]
################################################################################





################################################################################
def interlockcheck(intvl1,intvl2,pkover):
    #test if the two intvls are interlocked, each intvl: [l1,l2,r1,r2]
    #default pkover=0.2, max overlap of internal arms of interlocked DGs.
    i11,i12,i13,i14 = intvl1[0],intvl1[1],intvl1[2],intvl1[3]
    i21,i22,i23,i24 = intvl2[0],intvl2[1],intvl2[2],intvl2[3]
    if i21>=i11*pkover+i12*(1-pkover) and i22<=i14*pkover+i13*(1-pkover) and \
       i23>=i13*pkover+i14*(1-pkover) or i11>=i21*pkover+i22*(1-pkover) and \
       i12<=i24*pkover+i23*(1-pkover) and i13>=i23*pkover+i24*(1-pkover):
        return 1
    else: return 0
################################################################################





################################################################################
def mergegapmax(intervals, gapmax):
    #Merges a list of intervals with a gap smaller than a given cutoff.
    if not intervals: return []
    #intervals = [list(i) for i in intervals] #convert to list if tuple
    #intervals.sort(key=lambda x: x[0]) #sort by start only
    merged = [intervals[0]]
    for start1,end1 in intervals[1:]: #1: current; 0: last
        start0,end0 = merged[-1]
        if start1-end0 <= gapmax: merged[-1][1] = max(end0,end1)
        else: merged.append([start1,end1]) #Add new intvl
    return merged
################################################################################





################################################################################
def makesamples():
    #Focusing on 30 samples from SHARCLIP, ignoring PARIS, SHARC and eCLIP
    from itertools import product
    ENCORE,iPSC5 = ['HepG2','K562'],['iPSC','neuPC','neuron','astro','peri']
    RBPs = ['HNRNPC','POLR2A','RBFOX2','TDP43','TIA1'] 
    exps = ['assomul','indepmul']
    ENCORElist = [''.join(i) for i in product(ENCORE,RBPs,exps)] #ENCORE exps
    iPSClist = [''.join(i) for i in product(iPSC5,['HNRNPC'],exps)] #iPSC5 exps
    samples = ENCORElist+iPSClist; return samples
################################################################################






################################################################################
#common variables for direct import
cwd = '/Users/lu/Documents/lulab/projects/sharclip/_iPSsplice/'
all7lines = ['HepG2','K562','iPSC','neuPC','neuron','astro','peri']
iPS5lines = ['iPSC','neuPC','neuron','astro','peri']
samples = makesamples() #sample name, e.g. HepG2HNRNPCassomul
indices = {samples[i]:i+9 for i in range(len(samples))} #e.g. 
rmatsfolders = [
    "astro_vs_neuPC_outBAM_igvgtf_rmats", "astro_vs_neuron_outBAM_igvgtf_rmats",
    "astro_vs_peri_outBAM_igvgtf_rmats", "iPSC_vs_astro_outBAM_igvgtf_rmats",
    "iPSC_vs_neuPC_outBAM_igvgtf_rmats", "iPSC_vs_neuron_outBAM_igvgtf_rmats",
    "iPSC_vs_peri_outBAM_igvgtf_rmats", "neuPC_vs_neuron_outBAM_igvgtf_rmats",
    "neuPC_vs_peri_outBAM_igvgtf_rmats", "peri_vs_neuron_outBAM_igvgtf_rmats",
    "HepG2_vs_K562_outBAM_igvgtf_rmats"]
################################################################################

