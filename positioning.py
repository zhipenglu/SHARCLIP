"""
##### 0. Basic description. 
positioning.py. zhipengluchina@gmail.com, 2025-06-20. 
Compute midpoint, density and positioning of each arm for gap1 alignments for
short/long/all alignments (inter-midpoint distance cutoff=50). Export 12 files:
1. midpoints: all/short/long mid
2. density in 11nt windows: all/short/long density
3. positioning in 11 vs. 61nt windows, all/short/long positioning
4. midpoints for the contacts, instead of the arms

##### 1. Major functions:
1. struclib, CIGARdiv(CIGAR):
2. struclib, tointerval(line): these two to process BAM files
3. bam2midpoints(sam,sample,spanmax): export 4 midpoint bedgraphs and errorlines
4. mid2score(): returns density/positioning scores bedgraphs
5. a set of functions: makeMids, sortMids, runShort/Long/All/Contacts.
6. finally process the input file.
7. TODO: For gapm, divide to multiple segments.

##### 2. after the steps above, convert bedgraphs to bw to save space. 
cmd=/project/zhipengl_72/labshare/kentutils/bedGraphToBigWig
genome=/project/zhipengl_72/labshare/hg38/hg38_SHARCLIP_sizes.txt
for file in *positioning.bedgraph; \
do (salloc python removeOut.py $file &); done
for file in *positioning.bedgraph; \
do (salloc $cmd $file $genome ${file/bedgraph/bw} &); done

##### 3. Example command for a few lines.
script=/project/zhipengl_72/labshare/scripts/positioning.py 
script=/Users/lu/Dropbox/_scripts/positioning.py
python $script HepG2_HNRNPC_asso_gap1_PICALM.bam 50 1 11 61 &

##### 4. Actual runs on CARC. Make the mid, density and positioning bw files
/project/jianhuib_673/sharclip/sharclip_merged_gapped_reads/gap1
astro/HepG2/iPSC/K562/neuPC/neuron/peri_HNRNPC_asso/indep_gap1_sorted.bam
#use a slurm file to automatic processing.

##### 5. Running one at a time only to calculate density/positioning:
cd /project/jianhuib_673/sharclip/sharclip_merged_gapped_reads/gap1/positioning
salloc python /project/zhipengl_72/labshare/scripts/positioning.py \
../HepG2_TDP43_indep_gap1_sorted.bam 50 1 11 61 &

"""

################################################################################
##### 0. separately make a header fict for later use
def getgenome(sam):
    samf = pysam.AlignmentFile(sam, "rb")
    header = samf.header.to_dict()
    genome = {d["SN"]:d["LN"] for d in header["SQ"]}
    return genome
################################################################################





################################################################################
##### 1. export midpoints of gap1 arms to 4 bgs: all, short, long, and shortcx.
def bam2midpoints(sam,spanmax):
    #input: gap1 sorted bam; out bedgraph: chrom,mid,mid+1,count.
    #Export midpoints when pos<=l1 (left end) or new chrom to reduce memory use.    
    samf = pysam.AlignmentFile(sam, "rb"); errlines = []; p = sam[:-4]
    f1,f2 = open(p+'_allmid','w'),  open(p+'_shortmid','w')
    f3,f4 = open(p+'_longmid','w'), open(p+'_shortcxmid','w')
    m1,m2,m3 = {},{},{} #k=[chrom,mid],v=count; m1=m2+m3
    cx = {} #contact freq between neighbor peaks for short-dist gap1 alignments
    chroms,pos = [],0 #track current position
    header = samf.header.to_dict()
    for line in samf:
        line = line.tostring(); record = line.split(); chrom = record[2]
        try: intvl1,intvl2 = tuple(tointerval(line))
        except: errlines.append(line); continue
        l1,l2,r1,r2 = *intvl1[2:4],*intvl2[2:4]
        ml,mr = int(sum(intvl1[2:4])/2)-1,int(sum(intvl2[2:4])/2)-1 #left/right
        chroms,pos = [chrom] if not chroms else chroms, l1 if not pos else pos
        for k in [(chrom,ml),(chrom,mr)]:
            x = (chrom,int((ml+mr)/2))
            if mr-ml<= spanmax: m2[k] = m2.get(k,0)+1; cx[k] = cx.get(x,0)+1
            else: m3[k] = m3.get(k,0)+1
        for k in m2: m1[k] = m1.get(k,0) + m2[k] #merge m2 and m3 to m1.
        for k in m3: m1[k] = m1.get(k,0) + m3[k]
        
        r = range(0) #export a subset of the values, up to coord l1.
        if chrom != chroms[-1]:
            keys = list(m2.keys())[:-1] + list(m3.keys())[:-1]
            r = range(pos,max([k[1] for k in keys])+1)
            chroms.append(chrom); chrom = chroms[-2]; pos = l1 #old chr, new pos
        elif chrom == chroms[-1] and l1 > pos:  r = range(pos,l1); pos = l1

        for i in r: #export values behind pos
            k = (chrom,i); co = [chrom,str(i),str(i+1)]
            if k in m1: f1.write('\t'.join(co+[str(m1[k])])+'\n'); del m1[k]
            if k in m2: f2.write('\t'.join(co+[str(m2[k])])+'\n'); del m2[k]
            if k in m3: f3.write('\t'.join(co+[str(m3[k])])+'\n'); del m3[k]
            if k in cx: f4.write('\t'.join(co+[str(cx[k])])+'\n'); del cx[k]
    samf.close()
    for k in m1: f1.write('\t'.join([str(i) for i in [*k,k[1]+1,m1[k]]])+'\n') 
    for k in m2: f2.write('\t'.join([str(i) for i in [*k,k[1]+1,m2[k]]])+'\n')
    for k in m3: f3.write('\t'.join([str(i) for i in [*k,k[1]+1,m3[k]]])+'\n')
    for k in cx: f4.write('\t'.join([str(i) for i in [*k,k[1]+1,cx[k]]])+'\n')     
    f1.close(); f2.close(); f3.close(); f4.close()
    errf = open(sam[:-4]+'_errorlines.sam','w')
    errf.write(''.join(errlines)); errf.close() #error
    return header
################################################################################





################################################################################
##### 2. convert midpoints to density and positioning scores
def mid2score(mid,bedgd,bedgp,genome,mincov,winsize,peaksize):
    #input:  mid points bedgraph file. Needs to be sorted. 
    #output: bedgraphs of scores for density and positioning. Is output sorted?
    #use the score() function to calculate density and positioning scores
    midf = open(mid,'r') #midpoint bedgraph intermediate file
    bedgdf,bedgpf = open(bedgd,'w'),open(bedgp,'w') #density and positioning
    ms = [] #midpoints list [[chrom,m,count],...]
    chrom,m,c = '',0,0; prev = '' #previous chromosome
    while True:
        line = midf.readline(); x = line.split(); EOF = 1 if not line else 0
        bin1 = [0,0]; deck = deque()
        if EOF:
            if not ms: break # file is empty
            bin1 = [ms[0][1],ms[-1][1]]; deck = deque(ms)
        else:
            chrom,m,c = x[0],int(x[1]),int(x[3])
            if prev == '': prev = chrom; ms = [[chrom,m,c]]; continue # line 1
            if chrom != prev:
                bin1 = [ms[0][1],ms[-1][1]]; deck = deque(ms)
                chrom,prev = prev,chrom; ms = [[chrom,m,c]] #reset chrom.
            else: #same chromosome
                if not ms: ms.append([chrom,m,c]); continue
                else:
                    if m-ms[-1][1] <= winsize: ms.append([chrom,m,c]); continue
                    else: #desert encountered.
                        bin1 = [ms[0][1],ms[-1][1]]; deck = deque(ms)
                        ms = [[chrom,m,c]]
        
        #new method, compute two sets of densities first, then take ratios
        dscores,pscores = [[chrom,0,0,0]],[[chrom,0,0,0]] #density,positioning
        dpeak,dwin = {},{} #coverage in peak and window regions
        for (_,m,c) in deck: #(chrom,m,count)
            pl,pr = m-int((peaksize-1)/2),m+int((peaksize+1)/2)
            wl,wr = m-int((winsize-1)/2),m+int((winsize+1)/2)
            for i in range(max(pl,0),min(pr,genome[chrom])):
                dpeak[(chrom,i)] = dpeak.get((chrom,i),0)+c/peaksize
            for i in range(max(wl,0), min(wr,genome[chrom])):
                dwin[(chrom,i)] = dwin.get((chrom,i),0)+c/winsize
        for k in sorted(dpeak.keys()): #merge neighbors with identical values
            lastd,lastp = dscores[-1],pscores[-1]
            d,p = round(dpeak[k],3),round(dpeak[k]/dwin[k]*peaksize/winsize,3)
            if k[1] == lastd[2] and d == lastd[3]: dscores[-1][2]+=1
            else: dscores.append(list(k)+[k[1]+1,d])
            if k[1] == lastp[2] and p == lastp[3]: pscores[-1][2]+=1
            else: pscores.append(list(k)+[k[1]+1,p])
        for i in dscores[1:]: bedgdf.write('\t'.join([str(j) for j in i])+'\n')
        for i in pscores[1:]: bedgpf.write('\t'.join([str(j) for j in i])+'\n')
        if EOF: break
    midf.close(); bedgdf.close(); bedgpf.close()
    return
################################################################################





################################################################################
##### 3. all processing functions.
def makeMids(sam,spanmax):
    print(str(datetime.datetime.today())+"\tComputing mid points")
    header = bam2midpoints(sam,spanmax) #export 3 mid point bedgraphs
    print(str(datetime.datetime.today())+"\tFinished computing mid points")
    genome = {d["SN"]:d["LN"] for d in header["SQ"]}; return genome

def sortMids(sam):
    print(str(datetime.datetime.today())+"\tSorting midpoint bedgraphs")
    a = sam[:-4]+'_{}mid' #for [shortmid, longmid, allmid, shortcxmid]
    fs = a.format('short'),a.format('long'),a.format('all'),a.format('shortcx')
    for f in fs: 
        cmd = 'LC_COLLATE=c sort -k1,1 -k2,2n {0} >{1}'
        os.system(cmd.format(f,f+'.bedgraph')); os.system('rm {0}'.format(f))
    print(str(datetime.datetime.today())+"\tFinished sorting"); return 

def runShort(sam,genome,mincov,winsize,peaksize):
    print(str(datetime.datetime.today())+"\tAnalyzing short spans")
    a = sam[:-4]+'_long{}.bedgraph' #for shortmid,shortd,shortp 
    files = a.format('mid'),a.format('density'),a.format('positioning')
    mid2score(*files,genome,mincov,winsize,peaksize)
    print(str(datetime.datetime.today())+"\tFinished short spans"); return 

def runLong(sam,genome,mincov,winsize,peaksize):
    print(str(datetime.datetime.today())+"\tAnalyzing long spans")
    a = sam[:-4]+'_long{}.bedgraph' #for longmid,longd,longp 
    files = a.format('mid'),a.format('density'),a.format('positioning')
    mid2score(*files,genome,mincov,winsize,peaksize)
    print(str(datetime.datetime.today())+"\tFinished long spans"); return 

def runAll(sam,genome,mincov,winsize,peaksize): 
    print(str(datetime.datetime.today())+"\tAnalyzing all spans")
    a = sam[:-4]+'_all{}.bedgraph' #for allmid,alld,allp
    files = a.format('mid'),a.format('density'),a.format('positioning')
    mid2score(*files,genome,mincov,winsize,peaksize)
    print(str(datetime.datetime.today())+"\tFinished all spans"); return

def runContacts(sam,genome,mincov,winsize,peaksize):
    print(str(datetime.datetime.today())+"\tAnalyzing short contacts")
    a = sam[:-4]+'_shortcx{}.bedgraph' #for shortcx,shortcxd,shortcxp
    files = a.format('mid'),a.format('density'),a.format('positioning')
    mid2score(*files,genome,mincov,winsize,peaksize)
    print(str(datetime.datetime.today())+"\tFinished short contacts"); return 
################################################################################





################################################################################
##### 7. process the input file
if __name__ == '__main__':
    import sys, os, pysam, re, math, datetime, itertools
    from struclib import CIGARdiv, tointerval
    from collections import deque
    if len(sys.argv) < 2:
        sys.exit("python positioning.py sam spanmax mincov peaksize winsize")
    sam = sys.argv[1]
    spanmax = int(sys.argv[2])  #default 50, based on gap1 median spans.
    mincov = int(sys.argv[3])   #default 1, for positioning, not for density. 
    peaksize = int(sys.argv[4]) #default 11, for each arm. All types?
    winsize = int(sys.argv[5])  #default 61

    genome = getgenome(sam)
    makeMids(sam,spanmax) #print(genome)
    sortMids(sam)
    runContacts(sam,genome,mincov,winsize,peaksize)
    runShort(sam,genome,mincov,winsize,peaksize)
    runLong(sam,genome,mincov,winsize,peaksize)
    runAll(sam,genome,mincov,winsize,peaksize) #allmid to density/positioning
################################################################################








