"""
spliceloop.py, zhipengluchina@gmail.com, 2025-08-18
spliceloop.py merges *corr.py and *main.py to output figures and lists.
computes correlation DG vs. inclusion across cell lines.
Find DGs that contribute to skipping based on coverage and conservation.
Pay attention to these parameters in plotting. change depending on data scale. 
maxinc = 0.95; markerscale = 10; figlen = 10

Summary: major analysis for all mechanisms:
Test correlation among AS exons, or among cell types.
output lists, Manhattan plots, correlation volcanos, and heatmaps/arcs/ssmodels
################################################################################

1. AS, MATS.JC.txt, coords, inc levels. 
2. DGs, DGdiv.bb, info about location and counts
3. conservation, phastCons (phyloP) location and levels
4. structure models, from fasta, stability, etc. 

##### TODO: merge the two scripts, and export a long list.
in rewriting, use named dicts instead of hierarchical lists, tuples, and dicts

A consistent strategy for analyzing all mechanisms.
see spliceblockSE.py for other discussions on meta-analysis. 
The example plotting scripts are largely finished: 
bam2heat.py, heatmap view of contacts, in pairwise comparisons
gap2arc.py, arc view of contacts, using midpoints
DG2model.py, secondary structures of most abundant DGs

1. readrMATS(), outputs asmdict, splicing dict. Common first step. 
   key=('chrY','+',ES,EE,ESup,EEup,ESdn,EEdn)
   value={'astro':[1,0,585],'neuron':[0.998,0.0026,234]}, example
   Start with SE, rewrite to extend functionalities to other AS types. 

2. other steps using different methods.
   A. for loops, find DGs between introns vs. intra-intron, build models, cons
   B. for blockers, find DGs for SS5 and SS3E, build ss models, get phastcons
   C. for bridges, find cons regions, DGs near the bridges, build models, etc. 
   D. for switches, use intervaltree, find DGs for two exons, build models, etc 
   Starts with SE, extend functionalities to other AS types. 

3. correlation analysis can be done after these matching steps.
   A. for loops, expects negative correlation
   B. for blockers, also expects negative
   C. for bridges, may be either positive or negative, meaningful for single AS
   D. for switches, similar to blockers and loops.
   Starts with SE, also extend functionalities to other AS types.


Secondary structure models are made to guide mechanistic studies even though
they are not necessary for loop formation. Other related scripts: 
spliceloopcorr.py separates SE and DG cov info but does not use getphastcons().
spliceloopcorr.py: computes correlation DG vs. inclusion across cell lines. 
spliceloopmain.py reads phastCons, but not separate SE and DG cov values.
spliceloop_single.py tests correlation in each sample from one SE.MATS.JC.txt
spliceloop_single.py: plots scatter+LOWESS and spearman correlation among SE. 

Inter-intron DGs do not necessarily promote skipping. They may not block access
by the spliceosome. This is particularly interesting. What are their properties?
They may not be real duplexes, but rather protein-mediated proximity. We can
check splice site strength. See predictions from the DeltaSplice paper. 

PCCR_hg38.bed, contains 842258 records. 
Limitations of PCCR predictions:
1. only intronic from protein coding genes
2. constrained by phastCons
3. all-to-all pairs within 10kb range. ATE1 has an example exception, ~27kb

##### Writing the conclusions here: 
Using the dominance criteria, we identified ~100,000 DGs that likely contribute
to ~10,000 SEs in each cell type. These DGs are supported by data from ... cell
lines. ... of them connect deeply conserved regions, many of them span long
distances (histogram or ... distribution). We also validate a number of them by
ASOs and mutagenesis. To determine if they are dynamic across cell lines, we
measured them across 7 lineages, and here are examples of strong positives. 

##### Dominance criteria:
1. The first one or few DGs accounting for >50%/90% reads? May only need a few 
2. High in both asso and indep. A vague criterion but easy to set.  
3. For ones primarily in asso but not much in indep, RBP bridges?
4. For ones high in both, check structure stability.
5. Conservation, more likely to be critical 

##### input files:
IEI, not really necessary, can be rewritten ... 
DGdiv.bg bigbed
SE.MATS.JC.txt
phastCons bigwig
genome fasta


##### Basic procedure. Rewrite here. Are these outdated? Noted 2026-01-04
1. read IEI annotations, getIEI()
2. extract splicing data from all SE.MATS.JC.txt, getSEall()
3. identify major DGs connecting these IEI annotations. mainDGs()
4. present dominant DGs in one scatter plot, then examples with validations

Read all SE files together into one dict:
one IEI, each cell type, show inclusion level, and associated main DGs (1 file)
rank them, strength relative to top intra-intron DGs, local ones are likly to
have deeper coverage so it is not a fair comparison. Some of them have low
coverage, so there is too much background. Rank them by inclusion levels
This plot will be dense, so we will show some IEI examples and heatmaps. 

##### Criteria: 
1. inter-intron DGs, list all DGs within the coordinates 
2. either one or multiple, filter by number, take few with high coverage.
3. compare asso vs. indep and calculate structures.  
3. is there association with skipping levels across cell types? probably many
4. are they associated with diseases? not very likely but worth testing. 
5. are they associated with conservation? Useful to filter them. 

##### export the following fields, old design? 
1. skipped exon coordinates, chrom:start-end for rapid display on IGV
2. inclusion level, using our own calculation
3. DGs, one per line, for an SE, including l1,l2,r1,r2,cov and structure model
4. include DGID for comparison to DGs from other cell lines.

##### Example prep command:
cd /Users/lu/Desktop/sharclip; mamba activate crssant
script=~/Dropbox/_scripts/spliceloop.py
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_sorted_Birch
DGdiv=~/Desktop/sharclip/${pre}_T20_sorted_DGdiv.bb
phastbw=~/Documents/lulab/projects/hg38/hg38.phastCons100way.bw

##### tests on three input datasets: 
fa=~/igv/genomes/hg38/hg38.fa
fa=~/Documents/lulab/projects/hg38/fasta/chr4.fa
python $script $fa $DGdiv $phastbw loopDGs_vs_IEI_7celllines_ADD1
python $script $fa $DGdiv $phastbw loopDGs_vs_IEI_7celllines
python $script $fa $DGdiv $phastbw loopDGs_vs_IEI_7celllines_first1000

# use ghostscript (gs) to quickly rasterize large vector pdf
gs -dNOPAUSE -dBATCH -sDEVICE=pdfimage24 -r300 -o rasterized.pdf input.pdf
"""

#if k != ('chr4','+',2926641,2926675): continue #ADD1 E15.
#step1: get SE data using readrMATS(): 
#step2: get all DGs using which function? 

"""
need a complex data structure,
key = AScoord: (chrom,strand, 6-numbers)
value = DGs, each DGID binds a dict of coords, cov, etc. 

"""



################################################################################
##### 1. get coverage of each DG and each sample from DGdiv.bb
def getDGsep(asmdict,DGdiv,samples,indices):
    import pyBigWig
    #asmdict: key=(chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn)
    #value={cell:[incm,incs,sum_IJC_SJC],...}
    #convert to :{(chrom,strand,ES,EE):[bed12,incdict,DG1dict,DG2dict]}
    #incdict: {cell:inclusion_level}
    #DGdict: {(ID,coords):{sample:cov,...}}
    #DG1: intron1/2, DG2: inter-intron; DG1 is not used here. 
    #used a 5nt adjustment to enforce coverage in the introns, not exons. 
    def divcomma(string): return [int(i) for i in string.split(',')]
    def overlap(s1,e1,s2,e2): return s1<=e2 and s2<=e1
    f = pyBigWig.open(DGdiv,'r'); genome = f.chroms()
    IEIdict1 = {}
    for SE in IEIdict:
        chrom,strand,ES,EE = SE; record = IEIdict[SE][0].split()
        start,end = int(record[1]),int(record[2])
        if chrom not in genome or ES == start or EE == end: continue
        entries = f.entries(chrom,start,end)
        if not entries: continue
        DG1dict,DG2dict = {},{}
        for entry in entries: #each entry is one DG in bed12+30 format. 
            s1,e1,string = entry; data = string.split(); #print(data);
            ID,sizes,starts = data[0],divcomma(data[7]),divcomma(data[8])
            coords = (s1,s1+sizes[0],e1-sizes[1],e1)
            sampledict = {s:int(data[indices[s]]) for s in samples} #{s:cov}
            if overlap(start+5,ES-5,s1,s1+sizes[0]) or \
               overlap(start+5,ES-5,e1-sizes[1],e1) or \
               overlap(EE+5,end-5,s1,s1+sizes[0]) or \
               overlap(EE+5,end-5,e1-sizes[1],e1): #for intron1/2
                DG1dict[(ID,coords)] = sampledict
            if start+5<s1+sizes[0] and s1<ES-5 and \
               EE+5<e1 and e1-sizes[1]<end-5: #between introns
                DG2dict[(ID,coords)] = sampledict
        if DG2dict: IEIdict1[SE] = IEIdict[SE][0:2]+[DG1dict,DG2dict]
    return IEIdict1
################################################################################





################################################################################
##### 3. plot a scatter of DGs vs. IEI annotations, ranked by inclusion.
#plot a volcano for DGs that correlate with SE.
#HepG2 and K562 are outliers for ADD1. Focus on the most abundant ones. Why?
def plotcorr(IEIdict,exp,figname,volcanofig,volcanotxt):
    #IEIdict:{(chrom,strand,ES,EE):[bed12,incdict,DG1dict,DG2dict]}
    #incdict: {cell:inc}, e.g. {'iPSC':0.39,...} #1-7 cell types. no data: -1
    #DGdict: {(ID,coords):{sample:cov}} #all 30 samples
    #each SE needs at least 10 reads of IJC+SJC to be kept.
    post = "HNRNPC{}mul".format(exp) #postfix of the sample
    corrdict = {} #DGID:[SE,coords,inclist,covlist,c,-np.log10(p)]
    for SE in IEIdict:
        #1. for each IEI, get cells with inc levels.
        #2. get all DGs and only cells that have the inclusion levels
        #3. then plot all DGs for one SE event.
        #4. currently using only inter-intron DG2dict for normalization #

        #A. normalize against total for all DGs in one cell type in all samples.
        #if SE != ('chr4','+',2926641,2926675): continue #ADD1 SE. 
        bed12,incdict,DG1dict,DG2dict = IEIdict[SE]; print(DG2dict)
        covintrons = sum([int(i) for i in bed12.split()[10].split(',')])
        cells = [cell for cell in incdict if incdict[cell]>=0]
        totals = {c:sum([DG2dict[DG][c+post] for DG in DG2dict]) for c in cells}
        for DG in DG2dict:
            for c in cells:
                s = c+post
                if totals[c] >= 0.1: DG2dict[DG][s] = DG2dict[DG][s]/totals[c]
                else: DG2dict[DG][s] = 0

        #B. make a list of inclist and covlist for correlation testing. 
        inclist = [incdict[cell] for cell in cells]
        if all([inclist[0] == x for x in inclist]): continue #constant list
        for DG in DG2dict:
            if not DG2dict[DG]: continue
            covlist = [DG2dict[DG][cell+post] for cell in cells]
            if all([covlist[0] == x for x in covlist]): continue #constant list
            c,p = spearmanr(inclist,covlist) #calculate correlation
            corrdict[DG] = [SE,inclist,covlist,c,-np.log10(p)]
            
        """
        #C. plot each DG across cell types in a scatter connected by lines. 
        cellsort = [i[0] for i in sorted(incdict.items(),key=lambda x:x[1])]
        inclist = [incdict[cell] for cell in cellsort]
        if all([inclist[0] == x for x in inclist]): continue #constant list
        samplesx = [cell+"HNRNPC{}mul".format(exp) for cell in cellsort]
        fig = plt.figure(figsize=(2,2)) #print("DG2dict size:", len(DG2dict))
        #print(cellsort); print(inclist); print(DG2dict)
        for ID in DG2dict:
            covlist = [DG2dict[ID][s] for s in samplesx] #print(covlist)
            if all([covlist[0] == x for x in covlist]): continue #constant list
            c,p = pearsonr(inclist,covlist) #calculate pearson correlation
            if p > 0.05: continue #only plot positive/negative corr
            params = {'marker':'.','mfc':'r','mec':'r','ms':1,'ls':'--','lw':1}
            plt.plot(inclist,covlist,**params) 
        plt.xlim(0,1); plt.savefig(figname); plt.close()
        """

    #D. plot volcano for correlation and -log10(pvalues). Plot histograms later. 
    cs = [v[3] for v in corrdict.values()]; #print(cs) #if v[3]>=np.log10(20)
    ps = [v[4] for v in corrdict.values()]
    plt.scatter(cs,ps); plt.ylim(-0.2,5); plt.savefig(volcanofig); plt.close()
    print("Total DGs analyzed:", len(ps))
    print("DGs with p-value<=0.05:", len([p for p in ps if p>=np.log10(20)]))
    
    #E. export correlation data to a file:
    #[DGID,SE,coords,inclist,covlist,c,-np.log10(p),max(covlist)]
    txtf = open(volcanotxt,'w')
    for K in corrdict:
        SE,inclist,covlist,c,logp = corrdict[K]
        maxcov = str(round(max(covlist),5))
        #if logp <= np.log10(20): continue
        SE = ','.join([str(i) for i in SE])
        coords = ','.join([str(i) for i in K[1]])
        inclist = ','.join([str(round(i,5)) for i in inclist])
        covlist = ','.join([str(round(i,5)) for i in covlist])
        c,logp = str(round(c,5)),str(round(logp,5))
        l = [K[0],SE,coords,inclist,covlist,c,logp,maxcov]
        txtf.write('\t'.join(l)+'\n')
    txtf.close()
    return
################################################################################



################################################################################
##### 3. Manhattan scatter plot of inter-intronic DGs vs. IEIs, ranked by inc
#need to export all main DGs for each AS event, and other useful information. 
def exportmainDGs(IEIdict,fadict,minj,figname,figlen,markerscale):
    #IEIdict:{(chrom,strand,ES,EE):[bed12,SEinfo,DG1info,DG2info,[]]}
    #SEinfo: {cell:[name,IJC,SJC]}
    #DGinfo: {DGID:[l1,l2,r1,r2,cov]} #not separating samples here.
    #all SHARCLIP datasets were summed, including 30 samples.
    #each SE needs at least 10 reads to be kept.
    #each list of DGs for an SE should be >0 reads to be kept
    #DGlist: relative to total cov of DG2info, and max cov of DG1info.

    #A. get all the SElist and DGlist for each SE, calculate relative DG cov. 
    IEIsorted = [] #[[[SElist],[DGlist]],...]
    for SE in IEIdict:
        SEinfo,DG1info,DG2info = IEIdict[SE][1:4];
        SElist,DGlist = [],[[],[]]
        if SEinfo: 
            for cell in SEinfo:
                IJC,SJC = SEinfo[cell][1:3];
                if IJC+SJC>=minj: SElist.append(IJC/(IJC+SJC))
        if DG1info:
            maxDG = max([i[4] for i in DG1info.values()])
            total = sum([i[4] for i in DG2info.values()])
            for DGID in DG2info: #inter-intronic
                DGlist[0].append(DG2info[DGID][4]/total) #relative to sum(DG2)
                DGlist[1].append(DG2info[DGID][4]/maxDG) #relative to max(DG1)
        if SElist and DGlist:
            #DGlist = get80pct(DGlist,0.1) #Do not filter for now.
            IEIsorted.append([SElist,DGlist])
    print("IEIs passing filter for SE vs DG comparison:", len(IEIsorted))

    #B. for major loops, model the structures
    bedout = []
    for SE in IEIdict:
        chrom,strand,_,_ = SE
        print(SE)
        DG2info = IEIdict[SE][3] #inter-intronic DGs {DGID:[l1,l2,r1,r2,cov]}
        DG2sorted = sorted(DG2info.values(),key=lambda x:x[4],reverse=True)
        print("Models of top ranked DGs:")
        for DG in DG2sorted[:10]:
            seq,struc,mfe = ssmodel(chrom,DG[:4],strand,fadict)
            pairs = db2pairs(struc); print(chrom,DG[:4],seq,struc,mfe)
            sep = 0 #find location of &
            for i in range(len(struc)):
                if struc[i] == '&': sep = i; break
            for pair in sorted(pairs):
                pair[0] = pair[0]+DG[0] if pair[0]<sep else pair[0]+DG[2]
                pair[1] = pair[1]+DG[0] if pair[1]<sep else pair[1]+DG[2]
                bedout.append([chrom,str(pair[0]),str(pair[1])])
    string = '\n'.join(['\t'.join(i) for i in bedout])+'\n'
    f = open("ADD1_arcs.bed",'w')
    f.write('track graphType=arc\n'+string); f.close()
        
    """ example output: 
    chr4 [2926562, 2926592, 2926989, 2927015]
    TCTGTAACCTGATGGCTGTGACTGAATGCA&GATGTAAGTGCAGCCTCGGTTCAGAC
    ((((.((((.((.((((((.(((...((((&..))))))))))))))))))))))). -22.5
    """      
    
    #C. plot the fraction of coverage for each DG in each IEI.
    IEIsorted.sort()
    xs,ys,cs,ss = [],[],[],[] #xs:location; ys:frac_inter; ss:frac_intra
    for i in range(len(IEIsorted)):
        y = IEIsorted[i][1][0] #relative to total of DG2info
        z = IEIsorted[i][1][1] #relative to max of DG1info
        xs.extend([i]*len(y)); ys.extend(y);
        cs.extend([1-j for j in y]); ss.extend(z)
    print("Number of quantified DGs:", len(ys))
    print("Number of DGs with frac>=0.01:", len([i for i in ys if i>=0.01]))
    print("Number of DGs with frac>=0.05:", len([i for i in ys if i>=0.05]))
    fig = plt.figure(figsize=(figlen,1))
    plt.scatter(xs,ys,c=cs,cmap='gray',s=[j*markerscale for j in ss])
    plt.ylim(0,1); plt.savefig(figname); plt.close()    
    return
################################################################################





################################################################################
##### 4. filtering. no need now, given the color coding of relatve cov.
def get80pct(l1,mincov): 
    #l1: [[v1,v2],...].
    #v1: relative to total cov of DG2info; v2: relative to max DG2info
    #mincov is set to 0.1 by default, to remove low cov DGs.
    #Exception: if last one evaluated is <0.1, still keep it. 
    l1 = sorted(l1,reverse=True)
    total = sum([i[0] for i in l1]); target = 0.8*total; cumul = 0; l2 = []
    for i in l1:
        if i[0] < mincov: break
        cumul += i[0]; l2.append(i)
        if cumul > target: break
    if not l2: l2 = l1[0]
    return l2
################################################################################





################################################################################
##### 5. get average phastCons conservation of two DG arms for each IEI.
def getphastcons(IEIdict,phastbw):
    #after getDGall(), add conservation information to each DG. 
    #input: {(chrom,strand,ES,EE):[bed12,SEinfo,DG1info,DG2info,[]]}
    #DGinfo: {DGID:[l1,l2,r1,r2,cov]}]}, DG1: intron1/2, DG2: inter-intron
    #for each inter-intronic DG, calculate the average cons from phastcons.
    #plot the IEIs, and for each IEI, the distribution of average ...
    def divcomma(string): return [int(i) for i in string.split(',')]
    def overlap(s1,e1,s2,e2): return s1<=e2 and s2<=e1
    f = pyBigWig.open(phastbw,'r'); genome = f.chroms()
    for SE in IEIdict:
        chrom,strand,ES,EE = SE; bed12,SEinfo,DG1info,DG2info,_ = IEIdict[SE]
        if not DG2info: continue
        for ID in DG2info:
            l1,l2,r1,r2 = DG2info[ID][:4]
            wl = f.stats(chrom,l1,l2,exact=1)[0]; wl = 0 if wl is None else wl
            wr = f.stats(chrom,r1,r2,exact=1)[0]; wr = 0 if wr is None else wr
            cons = (wl*(l2-l1)+wr*(r2-r1))/(l2-l1+r2-r1)
            IEIdict[SE][4].append(cons)
    return IEIdict
################################################################################





################################################################################
##### 6. Manhattan plot of PhastCons 100-way scores for DGs in each IEI, with
#min inclusion scores <=0.95 for at least one sample.
def plotcons(IEIdict,minj,figname,figlen,markerscale):
    #IEIdict:{(chrom,strand,ES,EE):[bed12,SEinfo,DG1info,DG2info,conslist]}
    IEIsorted = [] #[[SElist,conslist],...]
    for SE in IEIdict:
        _,SEinfo,_,_,conslist = IEIdict[SE]; SElist = []#print(len(IEIdict[SE]))
        if SEinfo: 
            for cell in SEinfo:
                IJC,SJC = SEinfo[cell][1:3];
                if IJC+SJC>=minj: SElist.append(IJC/(IJC+SJC))
        if SElist and conslist: IEIsorted.append([SElist,conslist])
    print("IEIs passing filter for SE vs PhastCons comparison:", len(IEIsorted))
    IEIsorted.sort()
    xs,ys,cs,ss = [],[],[],[] #xs:location; ys:conslist; cs: colors; ss: sizes
    for i in range(len(IEIsorted)):
        y = IEIsorted[i][1]
        xs.extend([i]*len(y)); ys.extend(y); cs.extend([1-j for j in y])
    print("Number of quantified DGs:", len(ys))
    print("Number of DGs with cons>=0.5:", len([i for i in ys if i>=0.5]))
    print("Number of DGs with cons>=0.9:", len([i for i in ys if i>=0.9]))
    fig = plt.figure(figsize=(figlen,1))
    plt.scatter(xs,ys,c=cs,cmap='gray',s=0.5,lw=0)
    plt.ylim(0,1); plt.savefig(figname); plt.close()
    return 
################################################################################





################################################################################
##### 5. process files. 
if __name__ == "__main__":
    import os; os.environ['OPENBLAS_NUM_THREADS'] = '1' #override a bug on CARC
    import sys,pyBigWig
    from datetime import datetime
    from matplotlib import pyplot as plt
    from struclib import makesamples,readrMATS,readfa,ssmodel,db2pairs
    from struclib import cwd,all7lines,iPS5lines,rmatsfolders
    if len(sys.argv) < 5:
        sys.exit("python spliceloop.py fasta DGdiv phastbw prefix")
    fasta,DGdiv,phastbw,figname = sys.argv[1:6]
    figname1 = figname+"_main.pdf"
    figname2 = figname+"_cons.pdf" #phastCons scores.
    minj = 10; maxinc = 0.95; markerscale = 10 ; figlen = 10 #transcriptome. 
    #minj = 10; maxinc = 1; markerscale = 10; figlen = 2     #one gene
    rmatsfolders = ["astro_vs_neuron_outBAM_igvgtf_rmats"]   #test 

    #A. get genome fasta file
    print(str(datetime.today())+'\tExtracting fasta sequences')
    fadict = readfa(fasta); print("fasta chroms:", list(fadict.keys()))

    #B. get the splicing data.
    print(str(datetime.today())+'\tExtracting splicing data')
    asmdict = readrMATS(cwd,rmatsfolders,'SE',maxinc);
    print('Splicing dict:', len(asmdict))
    #5 cell lines, n=186102 for SE. After filtering with maxinc=0.98, n=40005.
    #for k in asmdict: print(k, asmdict[k]); break
    #key=('chrY','+',ES,EE,ESup,EEup,ESdn,EEdn)
    #value={'astro':[1,0,585],'neuron':[0.998,0.0026,234]}
    
    #C. get all the DG counts for each sample.
    print(str(datetime.today())+'\tExtracting DGs')
    IEIdict = getDGsep(asmdict,DGdiv); d = IEIdict
    print("Total number of inter-intron gap1 reads:",
          sum([0 if not d[SE][2] else sum([i[-1] for i in d[SE][2].values()])
               for SE in d]))

    #D. get all the PhastCons 100way data.
    print(str(datetime.today())+'\tExtracting PhastCons')
    IEIdict = getphastcons(IEIdict,phastbw); d = IEIdict
    
    #E. plot DGs for each IEI, sorted by inclusion levels.
    print(str(datetime.today())+'\tPlotting IEI DGs vs. splicing data (SE)')
    exportmainDGs(IEIdict,fadict,minj,figname1,figlen,markerscale)
    plotcons(IEIdict,minj,figname2,figlen,markerscale)

    #F. export all the AS, DG, structure, and conservation data to a table
    
    
    """
    #build models for all DGs? lets do it now    
    Total number of IEIs annotated: 194944
    IEIs with IJC+SJC>=10 and inc<=0.95 in >=1 samples 27149
    Number of quantified DGs: 1535085
    """
################################################################################




"""
No need to rerun this plotting, but will export all useful data. 

#all genes, took ~ 9 mins. 
python $script $IEI $DGdiv $phastbw loopDGs_vs_IEI_7celllines      
2025-08-23 13:49:42.314602	Extracting IEI annotations
Total number of IEIs annotated: 194944
2025-08-23 13:49:42.800200	Extracting splicing data
IEIs with IJC+SJC>=10 & inc<=0.95 in >=1 samples 27149
2025-08-23 13:49:45.235159	Extracting DGs
Total number of inter-intron gap1 reads: 554253693
2025-08-23 13:55:54.835175	Extracting PhastCons
2025-08-23 13:56:41.974163	Plotting IEI DGs vs. splicing data (SE)
Number of quantified DGs: 1535085
Number of DGs with frac>=0.01: 360821
Number of DGs with frac>=0.05: 77352
Number of quantified DGs: 1535085
Number of DGs with cons>=0.5: 82128
Number of DGs with cons>=0.9: 12745
"""

