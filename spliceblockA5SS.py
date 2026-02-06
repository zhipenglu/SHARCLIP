"""
spliceblockA5SS.py, 2025-10-02, zhipengluchina@gmail.com
Rewritten on 2025-10-04 from spliceblockSE.py. To be tested carefully. 
see parameters within script. Replacing previous two versions.

##### Basic procedure. 
1. read rMATS output AS data from the 5 cell lines, readrMATS() 
2. read DG from DGdiv.bb
3. build bp models for top DGs that block core motifs
4. extract conservation from phastCons
5. compute correlation among inc levels, gap1 reads, and structure stability. 

##### Major function:
A5SSblockers(asmdict,fadict,DGdiv,phastCons,exp,indices)
A. Get DGs in two regions: SS5 and BP-PPT-SS3
B. To compare genes in one sample, use MFE, no need for normalization? TO DO.
C. To compare one AS event across samples, use cov when MFE is OK. 

#conditions for blocking SS5:
A. For DGs, overlaps either arm, or in the loop and loop size<=50nts
B. For duplex models: overlaps either arm, or in the loop and loop size<=50nts

#conditions for blocking SS3 extended region (SS3 + 40nts, aka SS3X):
A. For DGs, overlaps either arm
B. For duplex models: overlaps either arm. What about loop region? Modify later.

visual inspection, unlikely to show strong correlation in MFE vs. inc levels.
should we expect strong association between structure strength and inc levels?
yes, but it depends on how the analysis is done.

##### Example A5SS from visual inspection, EEF1AKMT2 chr10:124788899-124789273

##### Example command. 
fasta=~/igv/genomes/hg38/hg38.fa
fasta=~/Documents/lulab/projects/hg38/fasta/chr16.fa
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs
DGdiv=${pre}_sorted_Birch_T20_sorted_DGdiv.bb
phastCons=~/Documents/lulab/projects/hg38/conservation/hg38.phastCons100way.bw
suffix=_all7lines_A5SSblockers_all.txt
python ~/Dropbox/_scripts/spliceblockA5SS.py $fasta $DGdiv $phastCons $suffix &



"""



################################################################################
##### 1. read DGs, find blockers, build models, compute correlation
def A5SSblockers(asmdict,fadict,DGdiv,phastCons,samples,
                 indices,flank,loopmax,DGext,topDGs,stdmin):
    """
    one function to get DGs, build models, and extract conservation.
    add new info dict to asmdict step by step.
    asmdict from readrMATS() in struclib.py
    input: asmdict[AScoord]['ASdict'] = {cell:[incm,incs,sum_IJC_SJC],...}
    output update:
    asmdict[AScoord]['covdicti'] = covdicti
    asmdict[AScoord]['duplexes3X'] = duplexes3X   #not needed here. 
    asmdict[AScoord]['duplexes5'] = duplexes5
    covdicti = {sample:count} #here count is the sum of HNRNPC asso and indep
    duplexes3X/5 = [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,
                   seq,struc,mfe,pairsg,label1,label2,cons1,cons2}

    loopmax=50, blocked in a tight loop, instead of the two arms
    flank=40, extending from SS3 to include PPT and BP. 
    DGext=10, extending each arm on each side by 10nts or until two arms merge.
    topDGs=10, top DGs by cov, the most likely blockers. 
    stdmin=0.1, inc level std min, 0.1 is decent variation based on prior tests
    """
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    def divcomma(string): return [int(i) for i in string.split(',')]
    DGdivf = pyBigWig.open(DGdiv,'r'); phastf = pyBigWig.open(phastCons,'r')
    for AScoord in asmdict: #process each AS event one by one. 
        chrom,strand,ESlong,EElong,ESshort,EEshort,ESflank,EEflank = AScoord
        if chrom not in DGdivf.chroms(): continue
        ASdict = asmdict[AScoord]['ASdict']
        std = np.std([i[0] for i in ASdict.values()])
        if std <= stdmin: continue #take most dramatic ones for initial tests. 



        #find example for A5SS for quick testing. #PKD1
        #if chrom != 'chr10': continue
        #if chrom != 'chr16' or ESlong>2113300 or EElong<2112700:continue
        

        
        #A. get gap1 cov in one intron as norm standards.
        intron = (EElong,ESflank) if strand=='+' else (EEflank,ESlong)
        entriesi = DGdivf.entries(chrom,*intron) #intron
        if not entriesi: continue
        covdicti = {s:0 for s in samples} #cov dict for two introns combined
        for entry in entriesi:
            s1,e1,string = entry; data = string.split()
            sizes,starts = divcomma(data[7]),divcomma(data[8])
            DGl,DGr = [s1,s1+sizes[0]],[e1-sizes[1],e1]
            if overlap(*DGl,*intron) or overlap(*DGr,*intron):
                for s in samples: covdicti[s] += int(data[indices[s]])
        asmdict[AScoord]['covdicti'] = covdicti
        
        #B. get all DGs overlapping SS5, a single nt position
        SS5long =  (chrom,EElong-1,EElong) if strand=='+' else \
                    (chrom,ESlong,ESlong+1) #for the long isoform
        SS5short = (chrom,EEshort-1,EEshort) if strand=='+' else \
                    (chrom,ESshort,ESshort+1)
        entries5long = DGdivf.entries(*SS5long)
        entries5short = DGdivf.entries(*SS5short)
        entries5long = [] if not entries5long else entries5long
        entries5short = [] if not entries5short else entries5short
        DGs5,duplexes5 = [],[] #all duplexes blocking SS5,
        #8 elements per DG up to now: [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict]
        entries5 = entries5long+entries5short
        if entries5:             
            for entry in entries5: #DGs
                s1,e1,string = entry; data = string.split(); DGID = data[0]
                sizes,starts = divcomma(data[7]),divcomma(data[8])                
                if s1<SS5long[1]<=s1+sizes[0] or e1-sizes[1]<=SS5long[1]<e1 or \
                   s1+sizes[0]<SS5long[1]<e1-sizes[1] and \
                   e1-sizes[1]-(s1+sizes[0])<=loopmax or \
                   s1<SS5short[1]<=s1+sizes[0] or e1-sizes[1]<=SS5short[1]<e1 \
                   or s1+sizes[0]<SS5short[1]<e1-sizes[1] and \
                   e1-sizes[1]-(s1+sizes[0])<=loopmax: #overlaps SS5long/short
                    covdict = {k:data[indices[k]] for k in indices}
                    c1 = [int(data[indices[k]]) for k in indices if asso in k]
                    c2 = [int(data[indices[k]]) for k in indices if indep in k]
                    nums = [s1,s1+sizes[0],e1-sizes[1],e1]
                    DGs5.append([DGID]+nums+[sum(c1),sum(c2),covdict])
            
            #C. get models for top DGs that overlap SS5
            if DGs5:
                #12 elements per DG now, c1/c2:HNRNPCasso/indep
                #[DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,seq,struc,mfe,pairsg]
                topDGs5 = sorted(DGs5,key=lambda x:sum(x[5:7]))[-topDGs:]
                for DG in topDGs5:
                    x = DGext #extend by 10nts on each side:
                    DG = [DG[0]]+[DG[1]-x,DG[2]+x,DG[3]-x,DG[4]+x]+DG[5:8]
                    if DG[2]>DG[3]: DG[2] = DG[3] = int(sum(DG[2:4])/2)
                    ss = seq,struc,mfe = ssmodel(chrom,DG[1:5],strand,fadict)
                    if mfe==0: duplexes5.append(DG+list(ss)+[[[0,0]]]); continue
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
                    if ssl1<SS5long[1]<ssl2 or ssr1<SS5long[1]<ssr2 or \
                       ssl2<=SS5long[1]<=ssr1 and ssr1-ssl2<=loopmax or \
                       ssl1<SS5short[1]<ssl2 or ssr1<SS5short[1]<ssr2 or \
                       ssl2<=SS5short[1]<=ssr1 and ssr1-ssl2<=loopmax:
                        duplexes5.append(DG+list(ss)+[pairsg])
                
                #D. get conservation for top duplexes (not DGs) overlapping SS5
                if duplexes5:
                    #16 elements now per DG: [DGID,l1,l2,r1,r2,sumc1,sumc2,
                    #covdict,seq,struc,mfe,pairsg,label1,label2,cons1,cons2]
                    for i in range(len(duplexes5)):
                        duplex = duplexes5[i]; pairsg = duplex[-1]
                        if pairsg[0] == [0,0]:
                            duplexes5[i] = duplex+['NA','NA',0,0]; continue
                        ssl1,ssl2,ssr1,ssr2 = sorted(pairsg[0]+pairsg[-1])                        
                        cons1 = phastf.stats(chrom,ssl1,ssl2,exact=True)[0]
                        cons2 = phastf.stats(chrom,ssr1,ssr2,exact=True)[0]
                        labels = []; x1 = SS5long[1:3]; x2 = SS5short[1:3]
                        if overlap(ssl1,ssl2,*x1): labels.append('SS5long')
                        elif overlap(ssl1,ssl2,*x2): labels.append('SS5short')
                        elif overlap(ssl1,ssl2,*intron): labels.append('I')
                        else: labels.append('E') #or O: other, e.g. exon
                        if overlap(ssr1,ssr2,*x1): labels.append('SS5long')
                        elif overlap(ssr1,ssr2,*x2): labels.append('SS5short')
                        elif overlap(ssr1,ssr2,*intron): labels.append('I')
                        else: labels.append('E')
                        duplexes5[i] = duplex+labels+[cons1,cons2]
        asmdict[AScoord]['duplexes5'] = duplexes5
        #E-G. get DGs, models and conservation for SS3X, BP+PPT+SS3, Not used.
    DGdivf.close(); phastf.close(); return asmdict #worked well. 
################################################################################





################################################################################
def exportasm(asmdict,samples,outfile): #export to a txt file
    """
    consider a general function for all the output formats. Not done yet. 
    SE and RI each has two sets of duplexes, SS5 and SS3X. 
    A3SS and A5SS only has one set of duplexes each. May be complicated. 
    
    how many changes in this script, compared to spliceblockSE.py:
    1. header 'SE' or others
    2. keys = 'ASdict','covdicti','duplexes3X','duplexes5', choose 3X and/or 5
    3. ASdict,covdicti,duplexes5 = ...
    4. duplexes = duplexes3X+duplexes5
    5. blocksite = 'SE' or others
    
    asmdict[AScoord]: {'ASdict':ASdict,'duplexes3X':[],'duplexes5':[]}
    no duplexes3X needed for A5SS
    gap1cv1/2 do not provide sufficient info for sorting across AS events,
    due to the variable coverage among introns. 
    """
    f = open(outfile,'w'); nDGs = 0
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    header = ['DGID','DGcoord','asso,indep','gap1asso,gap1indep','gap1total',
              'gap1cv1','gap1cv2','gap1block','gap1introns','seq','struc',
              'mfe','basepairs','A5SS','blocksite','incstd','incdict',
              'phastCons1','phastCons2','phastConsmean1','phastConsmean2']
    f.write('\t'.join(header)+'\n')
    for AScoord in asmdict:
        if 'covdicti' not in asmdict[AScoord]: continue
        keys = 'ASdict','covdicti','duplexes5'
        ASdict,covdicti,duplexes5 =[asmdict[AScoord][k] for k in keys]
        incstd = str(np.std([ASdict[c][0] for c in ASdict]))
        AScoordstr = ','.join([str(i) for i in AScoord])
        #note: not every AS event has inclusion data from all cell lines
        ASdictstr = ';'.join(['{}={}'.format(c,ASdict[c]) for c in ASdict])
        covistr = ';'.join(['{}={}'.format(s,covdicti[s]) for s in covdicti])
        #covistr, cov dict intron string, 30 samples, all DGs in two introns
        exps = 'HNRNPCassomul,HNRNPCindepmul' #experiments

        duplexes = duplexes5
        for i in range(len(duplexes)):
            nDGs+=1; blocksite = 'SS5'
            duplex = duplexes[i]
            DGID,l1,l2,r1,r2,sumc1,sumc2,covdict = duplex[:8]
            seq,struc,mfe,pairsg = duplex[8:12]
            label1,label2,cons1,cons2 = duplex[12:16]
            DGcoordstr = ','.join([str(i) for i in [l1,l2,r1,r2]])
            sums = str(sumc1)+','+str(sumc2); sumall = str(sumc1+sumc2)
            covstr = ';'.join(['{}={}'.format(s,covdict[s]) for s in samples])
            #covstr: gap1 cov dict string, 30 samples. 
            pairsstr = ','.join(['{}-{}'.format(*p) for p in pairsg])
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
    import numpy as np
    from matplotlib import pyplot as plt
    from struclib import readrMATS,readfa,ssmodel,db2pairs,overlap,makesamples
    from struclib import cwd,rmatsfolders,samples,indices

    if len(sys.argv)<4:
        sys.exit("python spliceblockA5SS.py fasta DGdiv phastCons")
    #4 types of input: fasta,DGdiv,phastCons,rmatsfolders
    #rmatsfolders = ["astro_vs_neuron_outBAM_igvgtf_rmats"] #fast test of one
    fasta,DGdiv,phastCons = sys.argv[1:4]
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

    #D. read the splicing data from rMATS output. ~34s for 10 SE comparisons. 
    #5 cell lines, n=186102 for SE. After filtering with maxinc=0.98, n=40005.
    asmdict = readrMATS(cwd,rmatsfolders,'A5SS',mininc,maxinc);
    print(str(datetime.datetime.today())+'\tFinished making splicing dict')
    
    #E. read DGs, find blockers, build models, extract phastCons
    #new data are added to the asmdict. 
    data = asmdict,fadict,DGdiv,phastCons,samples,indices
    params = flank,loopmax,DGext,topDGs,stdmin
    asmdict = A5SSblockers(*data,*params); #print(asmdict)
    print(str(datetime.datetime.today())+'\tFinished processing all data')

    #F. export data, 21 columns. Examples are analyzed later.
    nDGs = exportasm(asmdict,samples,outtxt)
    print(str(datetime.datetime.today())+'\tFinished exporting blockers')
    #for k in asmdict: print(asmdict[k]['ASdict']); break; sys.exit()
    #key=('chrY','+',ES,EE,ESup,EEup,ESdn,EEdn)
    #'ASdict': {'astro':[1,0,585],'neuron':[0.998,0.0026,234]},...}
    #in addition: covdict1, duplexes3X and duplexes5 lists

    print(); filters = (topDGs,mininc,maxinc,stdmin)
    print("Filters: topDGs={},mininc={},maxinc={},stdmin={}".format(*filters))
    print("Other: flank={},loopmax={},DGext={}".format(flank,loopmax,DGext))
    print("Number of AS events passing filter: {}".format(len(asmdict)))
    print('Exported DG blockers: {}\n'.format(nDGs))
    
    #G. computes correlation and plot correlation scatters 
    #see details in a separate script: splicevsDG.py
################################################################################



"""
Are there any published examples of structures regulating A3SS or A5SS?
See Shepard and Hertel 2008. Three sup tables included a small number of
potentially conserved structures covering A3SS, A5SS, and SE.

Coordinates were in hg18. Lifted to hg38 as follows:
chrom Evofold_start Evofold_stop ss1 ss2 name
chr1 154153118 154153172 154153123 154153151 KIAA0907
The regulatory effects will be most obvious for alt. splice sites that are
distant from each other. For those that are too close the blocking will affect
both, likely to the same extent.

hg38, all 7 cell lines. 
2025-10-13 23:18:55.364370	Starting analysis
2025-10-13 23:19:12.763782	Finished reading fasta
2025-10-13 23:19:16.617270	Finished making splicing dict
2025-10-14 00:15:52.654393	Finished processing all data
2025-10-14 00:15:58.262787	Finished exporting blockers

Filters: topDGs=10,mininc=0,maxinc=1,stdmin=0
Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 15716
Exported DG blockers: 123762
"""


