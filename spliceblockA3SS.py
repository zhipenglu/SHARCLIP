"""
spliceblockA3SS.py, 2025-10-02, zhipengluchina@gmail.com
Rewritten on 2025-10-03 from spliceblockSE.py. To be tested carefully. 
See parameters within script.
Master output is done. Next test and find correlations. 

##### Major considerations:
1. Find DGs that cover the two A3SS regions, many of which may overlap.
2. See plot_inclusion.py for histogram of distance distributions.
3. Based on this consideration, we will treat the alt. sites separately.

##### Possibilities to consider in labeling the mechanisms (not separated):
1. blocking SS3 for both at the same time
2. blocking one SS3 and the BP+PPT for the other.
3. blocking BP+PPT for both
4. blocking SS3 or BP+PPT for one, and none for the other.
5. blocking SS3 or BP+PPT for one, and looping the other.
6. other possibilities exist. 

##### Changes from spliceblockSE.py:
1. only duplexes3, but not duplexes5. Also only duplexes5 for A5SS later. 

##### Basic procedure. 
1. read rMATS output AS data from the 5/7 cell lines, readrMATS() 
2. read DG from DGdiv.bb
3. build bp models for top DGs that block core motifs
4. extract conservation from phastCons
5. compute correlation among inc levels, gap1 reads, and structure stability. 

##### Comment on 2025-09-29. A major technical challenge in this study is the
complexity of data structures. We need to simultaneously process multiple input
and output data types, including fasta, DG coords from bigbed, DG coverage
dictionary, AS coords, AS inclusion levels in dicts, and conservation bigwig
(phastCons). Here I describe the framework for standardizing data structures. 

##### Major function:
A3SSblockers(asmdict,fadict,DGdiv,phastCons,exp,indices)
A. Get DGs in two regions: SS5 and BP-PPT-SS3
B. To compare genes in one sample, use MFE, no need for normalization? TO DO.
C. To compare one AS event across samples, use cov when MFE is OK. 
#conditions for blocking SS5:
A. For DGs, overlaps either arm, or in the loop and loop size<=50nts
B. For duplex models: overlaps either arm, or in the loop and loop size<=50nts
#conditions for blocking SS3 extended region (SS3 + 40nts, SS3X):
A. For DGs, overlaps either arm
B. For duplex models: overlaps either arm. What about loop region? Modify later.

Find a different gene for A3SS or A5SS. 
##### Example: focus on IDE for a quick test, for SE
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs
DGdiv=${pre}_sorted_Birch_T20_sorted_DGdiv.bb
phastCons=~/Documents/lulab/projects/hg38/conservation/hg38.phastCons100way.bw
fasta=~/igv/genomes/hg38/hg38.fa
fasta=~/Documents/lulab/projects/hg38/fasta/chrY.fa
suffix=_all7lines_A3SSblockers_all.txt
python ~/Dropbox/_scripts/spliceblockA3SS.py $fasta $DGdiv $phastCons $suffix &

"""



################################################################################
##### 1. read DGs, find blockers, build models, compute correlation
def A3SSblockers(asmdict,fadict,DGdiv,phastCons,samples,
                 indices,flank,loopmax,DGext,topDGs,stdmin):
    """
    one function to get DGs, build models, and extract conservation.
    add new info dict to asmdict step by step.
    asmdict from readrMATS(), from struclib.py 
    input: asmdict[AScoord]['ASdict'] = {cell:[incm,incs,sum_IJC_SJC],...}
    output update:
    asmdict[AScoord]['covdicti'] = covdicti
    asmdict[AScoord]['duplexes3X'] = duplexes3X
    asmdict[AScoord]['duplexes5'] = duplexes5 #not needed here
    covdicti = {sample:count} #here count is the sum of HNRNPC asso and indep
    duplexes3X/5 = [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,
                    seq,struc,mfe,pairsg,label1,label2,cons1,cons2} #modify

    loopmax=50, blocked in a tight loop, instead of the two arms
    flank=40, extending from SS3 to include PPT and BP. 
    DGext=10, extending each arm on each side by 10nts or until two arms merge.
    topDGs=10, top DGs by cov, the most likely blockers. 
    mininc=0; maxinc=1, parameters used in readrMATS(), not here. 
    stdmin=0, inc level std min, 0.1 is decent variation based on prior tests
    """
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    def divcomma(string): return [int(i) for i in string.split(',')]
    DGdivf = pyBigWig.open(DGdiv,'r'); phastf = pyBigWig.open(phastCons,'r')
    for AScoord in asmdict: #process each AS event one by one. 
        chrom,strand,ESlong,EElong,ESshort,EEshort,ESflank,EEflank = AScoord
        if chrom not in DGdivf.chroms(): continue
        ASdict = asmdict[AScoord]['ASdict']
        std = np.std([i[0] for i in ASdict.values()])
        if std <= stdmin: continue #take most dramatic ones for statistics.


        
        #it seems complete switching is very rare for A3SS/A5SS, why? ...        
        #if chrom != 'chr1' or ES>225520000 or EE<225514500: continue #ENAH
        #if chrom != 'chrY': continue


        
        #A. get gap1 cov in one intron as norm standards (min intron).
        intron = (EElong,ESflank) if strand=='-' else (EEflank,ESlong)
        print(chrom,*intron)
        entriesi = DGdivf.entries(chrom,*intron) #intron
        if not entriesi: continue
        covdicti = {s:0 for s in samples} #cov dict for two introns combined
        for entry in entriesi:
            s1,e1,string = entry; data = string.split()
            sizes,starts = divcomma(data[7]),divcomma(data[8])
            DGl,DGr = [s1,s1+sizes[0]],[e1-sizes[1],e1] #left or right arm
            if overlap(*DGl,*intron) or overlap(*DGr,*intron):
                for s in samples: covdicti[s] += int(data[indices[s]])
        asmdict[AScoord]['covdicti'] = covdicti
        #B-D. get DGs, models and conservation for SS5. Not needed here. 
        
        #E. get all DGs overlapping SS3X or BP+PPT+SS3, e.g. a 40nt region
        SS3Xlong =  (chrom,ESlong-flank,ESlong) if strand=='+' else \
                    (chrom,EElong,EElong+flank) #for the long isoform
        SS3Xshort = (chrom,ESshort-flank,ESshort) if strand=='+' else \
                    (chrom,EEshort,EEshort+flank)
        entries3Xlong = DGdivf.entries(*SS3Xlong)
        entries3Xshort = DGdivf.entries(*SS3Xshort)
        entries3Xlong = [] if not entries3Xlong else entries3Xlong
        entries3Xshort = [] if not entries3Xshort else entries3Xshort
        DGs3X,duplexes3X = [],[] #all duplexes blocking SS3.
        #12 elements per DG now: [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict]
        entries3X = set(entries3Xlong+entries3Xshort)
        if entries3X:
            for entry in entries3X: #DGs
                s1,e1,string = entry; data = string.split(); DGID = data[0];
                sizes,starts = divcomma(data[7]),divcomma(data[8])
                if overlap(s1,s1+sizes[0],*SS3Xlong[1:]) or \
                   overlap(s1,s1+sizes[0],*SS3Xshort[1:]) or \
                   overlap(e1-sizes[1],e1,*SS3Xlong[1:]) or \
                   overlap(e1-sizes[1],e1,*SS3Xshort[1:]): #SS3Xlong or short
                    covdict = {k:data[indices[k]] for k in indices}
                    c1 = [int(data[indices[k]]) for k in indices if asso in k]
                    c2 = [int(data[indices[k]]) for k in indices if indep in k]
                    nums = [s1,s1+sizes[0],e1-sizes[1],e1]
                    DGs3X.append([DGID]+nums+[sum(c1),sum(c2),covdict])
                        
            #F. get models for top DGs overlapping SS3X
            if DGs3X:
                #12 elements per DG now, c1/c2: HNRNPCasso/indep
                #[DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,seq,struc,mfe,pairsg]
                topDGs3X = sorted(DGs3X,key=lambda x:sum(x[5:7]))[-topDGs:]
                for DG in topDGs3X:
                    x = DGext; #extend by 10nts on each side:
                    DG = [DG[0]]+[DG[1]-x,DG[2]+x,DG[3]-x,DG[4]+x]+DG[5:8]
                    if DG[2]>DG[3]: DG[2] = DG[3] = int(sum(DG[2:4])/2)
                    ss = seq,struc,mfe = ssmodel(chrom,DG[1:5],strand,fadict)
                    if mfe==0: duplexes3X.append(DG+list(ss)+[[[0,0]]]);continue
                    pairs = db2pairs(struc); len1 = len(struc.split('&')[0])                    
                    ssl1,ssl2 = pairs[0][0]+DG[1],pairs[-1][0]+DG[1]
                    ssr1,ssr2 = pairs[-1][1]+DG[3]-len1,pairs[0][1]+DG[3]-len1
                    ssl1,ssl2,ssr1,ssr2 = sorted([ssl1,ssl2,ssr1,ssr2])
                    pairsg = [] #get base pairs in genomic coordinates
                    sep = list(struc).index('&')
                    for pair in sorted(pairs):
                        p0,p1 = pair; pair[0] = p0+DG[1] if p0<sep else p0+DG[3]
                        pair[1] = p1+DG[1]-len1 if p1<sep else p1+DG[3]-len1
                        pairsg.append(pair)
                    if overlap(ssl1,ssl2,*SS3Xlong[1:]) or \
                       overlap(ssr1,ssr2,*SS3Xlong[1:]) or \
                       overlap(ssl1,ssl2,*SS3Xshort[1:]) or \
                       overlap(ssr1,ssr2,*SS3Xshort[1:]):
                        duplexes3X.append(DG+list(ss)+[pairsg])

                #G. get conservation for top duplexes (not DG) that overlap SS3X
                if duplexes3X:
                    #16 elements per DG now, [DGID,l1,l2,r1,r2,sumc1,sumc2,
                    #covdict,seq,struc,mfe,pairsg,label1,label2,cons1,cons2]
                    for i in range(len(duplexes3X)):
                        duplex = duplexes3X[i]; pairsg = duplex[-1]
                        if pairsg[0]==[0,0]:
                            duplexes3X[i] = duplex+['NA','NA',0,0]; continue
                        ssl1,ssl2,ssr1,ssr2 = sorted(pairsg[0]+pairsg[-1])
                        cons1 = phastf.stats(chrom,ssl1,ssl2,exact=True)[0]
                        cons2 = phastf.stats(chrom,ssr1,ssr2,exact=True)[0]
                        labels = []; x1 = SS3Xlong[1:3]; x2 = SS3Xshort[1:3]
                        if overlap(ssl1,ssl2,*x1): labels.append('SS3Xlong')
                        elif overlap(ssl1,ssl2,*x2): labels.append('SS3Xshort')
                        elif overlap(ssl1,ssl2,*intron): labels.append('I')
                        else: labels.append('E') #or O: other, e.g. exon
                        if overlap(ssr1,ssr2,*x1): labels.append('SS3Xlong')
                        elif overlap(ssr1,ssr2,*x2): labels.append('SS3Xshort')
                        elif overlap(ssr1,ssr2,*intron): labels.append('I')
                        else: labels.append('E')
                        duplexes3X[i] = duplex+labels+[cons1,cons2]
        asmdict[AScoord]['duplexes3X'] = duplexes3X
    DGdivf.close(); phastf.close(); return asmdict #worked well. 
################################################################################





################################################################################
def exportasm(asmdict,samples,outfile): #export to a txt file
    """
    asmdict[AScoord]: {'ASdict':ASdict,'duplexes3X':[],'duplexes5':[]}
    no duplexes5 needed for A3SS. 
    ASdict: {cell:[incm,incs,sum_IJC_SJC],...}
    should I use the best structure or the integration of all top structures?
    the original structures of each DG is likely to be incomplete.
    for correlation among AS events, how do we normalize? see spliceloop.py
    output 21 columns in 5 groups: DGout, gap1out, ssout, ASout, consout

    gap1cv1/2 do not provide sufficient info for sorting across AS events,
    due to the variable coverage among introns. 
    """
    f = open(outfile,'w'); nDGs = 0
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    header = ['DGID','DGcoord','asso,indep','gap1asso,gap1indep','gap1total',
              'gap1cv1','gap1cv2','gap1block','gap1introns','seq','struc',
              'mfe','basepairs','A3SS','blocksite','incstd','incdict',
              'phastCons1','phastCons2','phastConsmean1','phastConsmean2']
    f.write('\t'.join(header)+'\n')
    for AScoord in asmdict:
        if 'covdicti' not in asmdict[AScoord]: continue
        keys = 'ASdict','covdicti','duplexes3X'
        ASdict,covdicti,duplexes3X = [asmdict[AScoord][k] for k in keys]
        incstd = str(np.std([ASdict[c][0] for c in ASdict]))
        AScoordstr = ','.join([str(i) for i in AScoord])
        #note: not every AS event has inclusion data from all cell lines
        ASdictstr = ';'.join(['{}={}'.format(c,ASdict[c]) for c in ASdict])
        covistr = ';'.join(['{}={}'.format(s,covdicti[s]) for s in covdicti])
        #covistr, cov dict intron string, 30 samples, all DGs in two introns
        exps = 'HNRNPCassomul,HNRNPCindepmul' #experiments

        duplexes = duplexes3X
        for i in range(len(duplexes)):
            nDGs+=1; blocksite = 'SS3X'; duplex = duplexes[i]
            DGID,l1,l2,r1,r2,sumc1,sumc2,covdict = duplex[:8]
            seq,struc,mfe,pairsg = duplex[8:12]
            label1,label2,cons1,cons2 = duplex[12:16]
            DGcoordstr = ','.join([str(i) for i in [l1,l2,r1,r2]])
            sums = str(sumc1)+','+str(sumc2); sumall = str(sumc1+sumc2)
            covstr = ';'.join(['{}={}'.format(s,covdict[s]) for s in samples])
            #covstr: gap1 cov dict string, 30 samples. 
            pairsstr = ','.join(['{}-{}'.format(*p) for p in pairsg])
            #for k in covdicti: print(covdict[k],covdicti[k]); 
            rel1 = [int(covdict[k])/covdicti[k] for k in covdicti if asso in k
                    and covdicti[k]] #relative coverage to the intron 
            cv1 = '0' if not sum(rel1) else str(np.std(rel1)/np.mean(rel1))
            rel2 = [int(covdict[k])/covdicti[k] for k in covdicti if indep in k
                    and covdicti[k]]
            cv2 = '0' if not sum(rel2) else str(np.std(rel2)/np.mean(rel2))
            DGout = [DGID,DGcoordstr] #n=2
            gap1out = [exps,sums,sumall,cv1,cv2,covstr,covistr] #n=7
            ssout = [seq,struc,str(mfe),pairsstr] #n=4
            ASout = [AScoordstr,blocksite,incstd,ASdictstr] #n=4
            consout = [label1,label2,str(cons1),str(cons2)] #n=4
            f.write('\t'.join(DGout+gap1out+ssout+ASout+consout)+'\n') #n=21
    f.close()
    return nDGs
################################################################################



    
    
################################################################################
##### process input files, merge several functions
if __name__ == '__main__':
    import sys, pyBigWig, datetime
    import numpy as np; from matplotlib import pyplot as plt
    from struclib import readrMATS,readfa,ssmodel,db2pairs,overlap,makesamples
    from struclib import cwd,rmatsfolders,samples,indices

    if len(sys.argv)<4:
        sys.exit("python spliceblockA3SS.py fasta DGdiv phastCons")
    #4 types of input: fasta,DGdiv,phastCons,rmatsfolders
    #rmatsfolders = ["astro_vs_neuron_outBAM_igvgtf_rmats"] #fast test of one
    fasta,DGdiv,phastCons,suffix = sys.argv[1:5]
    outtxt = DGdiv.split('.')[0]+suffix

    #A. parameters to filter DGs, AS inclusion and variation among cell types
    #vary the topDGs, mininc and maxinc levels to find more or less DGs
    #mininc, i.e., highest inc levels among cell lines should be >=mininc
    #maxinc, i.e., lowest inc levels among cell lines should be <=maxinc
    topDGs = 10 #int(1E9) #top DG blockers based on HNRNPC asso+indep coverage
    mininc,maxinc = 0,1
    stdmin = 0 #std of inc levels across cell lines. 

    #B. parameters not meant to change. 
    flank = 40 #range from SS3, as the SS3X (extended region)
    loopmax = 50 #max len of loops in DGs and duplexes that act as blockers
    DGext = 10 #extension of DG arms in 4 directions, without causing overlaps
    
    #C. read genome fasta file, ~16s for hg38.fa
    print(str(datetime.datetime.today())+'\tStarting analysis')
    fadict = readfa(fasta) 
    print(str(datetime.datetime.today())+'\tFinished reading fasta')

    #D. read the splicing data from rMATS output. ~...s for 10 A3SS comparisons. 
    #5 cell lines, n=186102 for SE. After filtering with maxinc=0.98, n=40005.
    #change these numbers for A3SS
    asmdict = readrMATS(cwd,rmatsfolders,'A3SS',mininc,maxinc);
    print(str(datetime.datetime.today())+'\tFinished making splicing dict')
    
    #E. read DGs, find blockers, build models, extract phastCons
    #new data are added to the asmdict. 
    data = asmdict,fadict,DGdiv,phastCons,samples,indices
    params = flank,loopmax,DGext,topDGs,stdmin
    asmdict = A3SSblockers(*data,*params); #print(asmdict)
    print(str(datetime.datetime.today())+'\tFinished processing all data')

    #F. export data, 21 columns. Examples are analyzed later.
    nDGs = exportasm(asmdict,samples,outtxt)
    print(str(datetime.datetime.today())+'\tFinished exporting blockers')
    #for k in asmdict: print(asmdict[k]['ASdict']); break; sys.exit()
    #key=('chrY','+',ES,EE,ESup,EEup,ESdn,EEdn). This example is SE. 
    #'ASdict': {'astro':[1,0,585],'neuron':[0.998,0.0026,234]},...}
    #in addition: covdict1, duplexes3X and duplexes5 lists

    print(); filters = (topDGs,mininc,maxinc,stdmin)
    print("Filters: topDGs={},mininc={},maxinc={},stdmin={}".format(*filters))
    print("Other: flank={},loopmax={},DGext={}".format(flank,loopmax,DGext))
    print("Number of AS events passing filter: {}".format(len(asmdict)))
    print('Exported DG blockers: {}\n'.format(nDGs))
    
    #G. computes correlation and plot correlation scatters 
    #see details in a separate script: spliceblockcorr.py
################################################################################



"""
two cell lines, hg38. 
2025-10-12 23:13:48.986713	Starting analysis
2025-10-12 23:14:05.307353	Finished reading fasta
2025-10-12 23:14:05.691178	Finished making splicing dict
2025-10-12 23:55:53.697506	Finished processing all data
2025-10-12 23:55:58.190019	Finished exporting blockers

Filters: topDGs=10,mininc=0,maxinc=1,stdmin=0
Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 13042
Exported DG blockers: 102949

All 7 cell lines. hg38
2025-10-13 13:17:24.164478	Starting analysis
2025-10-13 13:17:43.831212	Finished reading fasta
2025-10-13 13:17:48.296480	Finished making splicing dict
2025-10-13 15:21:45.878112	Finished processing all data
2025-10-13 15:21:52.391180	Finished exporting blockers

Filters: topDGs=10,mininc=0,maxinc=1,stdmin=0
Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 17818
Exported DG blockers: 138950

DONE? ANY PROBLEMS? 

"""


