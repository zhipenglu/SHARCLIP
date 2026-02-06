"""
spliceblockSE.py, 2025-10-02, zhipengluchina@gmail.com
see parameters within script. Replacing previous two versions.
Master output is done. Test and find correlations using splicevsDG_celllines.py
Next iteration:
topDGs to 1E9, effectively removing the limit. 
maxinc to 0.95, including more SE events. 
cell lines to 7, including all the datasets

Plots based on the complete blockers table: DONE.  
1. plot DG average PhastCons of the intronic one or both, or the intronic one.
2. plot fractions of reads in each DG? or MFE? or both?
3. plot scatter of MFE vs. fraction to visualize diff between asso and indep 
4. summary plots of total numbers, DGs per AS, etc. see plotmechs.py

##### Basic procedure. 
1. read rMATS output AS data from the 7 cell lines. readrMATS() from struclib
2. read DGs from DGdiv.bb
3. build bp models for top DGs that block core motifs
4. extract conservation from phastCons
5. export blockers to a txt file, 21 columns
6. splicevsDG_celllines.py: correlate PSI, gap1 reads, and structure stability. 

##### Comment on 2025-09-29. A major technical challenge in this study is the
complexity of data structures. We need to simultaneously process multiple input
and output data types, including fasta, DG coords from bigbed, DG coverage
dictionary, AS coords, AS inclusion levels in dicts, and conservation bigwig
(phastCons). Here I describe the framework for standardizing data structures. 
should I use the best structure or the integration of all top structures?
the original structures of each DG is likely to be incomplete.

For correlation among AS events, how do we normalize? For blockers, see
splicevsDG_celllines.py, neighbor exons for SE/MXE and A5SS. Tested the intron
for RI. The correlations for SE were very strong, but not for A5SS or RI. A3SS
is ignored due to the close proximity of BP/PPT/SS3 for both splice sites. 

##### General output considerations for all AS structural mechanisms:
1. Key information: DG, AS, structure, and conservation, across all samples
2. Other info, e.g. disease variants, function annotations are considered later. 
3. One DG per line for blockers, loops, and bridges, two per line for switches
4. switches are not analyzed systematically now. Only mannual curation for MXE.
4. sorting: DG (rel) cov total, cov std, MFE, AS inc std, arm/nearby cons, etc.
5. header and meta-info lines should be included for clarity and publication

##### List of columns, also applicable to loops and bridges, but not switches
1. DG ID, 4 elements, gene,gene,num,frac. 
2. DG coords, 4 numbers, list sep by comma
3. DG name of the multiple samples, e.g. HNRNPCassomul,HNRNPCindepmul
4. DG cov total asso and indep
5. DG cov total, HNRNPC asso+indep
6. DG cov coefficient of variation, after norm, for selecting dynamic DGs, asso
7. DG cov coefficient of variation, after norm, for selecting dynamic DGs, indep
8. DG cov dict, e.g. sample=value; ..., ~30 samples
9. DG cov dict, normalization reference, needed for stats analysis. 

10. struc, seq NNN&NNN, length depends on DG size, non-uniform, '+' strand
11. struc, dot-bracket format, (((&))), '+' strand
12. struc MFE, useful for quick sorting
13. struc, pairs of numbers, e.g. x1-y1,x2-y2, for conversion to bed arcs

14. AS coord, 8 elements, chrom,strand and 6 numbers for SE/A3SS/A5SS/RI.
15. AS blocker type, e.g. SS5 or SS3X (extended)
16. AS incstd, useful for quick sorting
17. AS inclusion dict, e.g. cell=value; ..., 5 or 7 cell lines

18. block location for left arm
19. block location for right arm
20. phastCons mean left arm of duplex
21. phastCons mean right arm of duplex

##### Major function:
SEblockers(asmdict,fadict,DGdiv,phastCons,exp,indices)
A. Get DGs in two regions: SS5 and BP-PPT-SS3
B. To compare genes in one sample, use MFE, no need for normalization? TO DO.
C. To compare one AS event across samples, use cov when MFE is OK. 
#conditions for blocking SS5:
A. For DGs, overlaps either arm, or in the loop and loop size<=50nts
B. For duplex models: overlaps either arm, or in the loop and loop size<=50nts
#conditions for blocking SS3 extended region (SS3 + 40nts, aka SS3X):
A. For DGs, overlaps either arm
B. For duplex models: overlaps either arm. What about loop region? Modify later.

##### Example: focus on TCF3 and IDE for quick tests
cd /Users/lu/Documents/lulab/projects/sharclip; mamba activate crssant
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs
DGdiv=${pre}_sorted_Birch_T20_sorted_DGdiv.bb
phastCons=~/Documents/lulab/projects/hg38/conservation/hg38.phastCons100way.bw
fasta=~/igv/genomes/hg38/hg38.fa
fasta=~/Documents/lulab/projects/hg38/fasta/chr21.fa
suffix=_all7lines_SEblockers_chr21.txt
script=~/Dropbox/_scripts/spliceblockSE.py
time python $script $fasta $DGdiv $phastCons $suffix &




##### Post processing, export base pairs. Example. Is it correct? 
f=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_sorted_Birch_\
T20_sorted_DGdiv_5iPSlines_SEblockers.txt
awk '$1~/7031/' $f | cut -f13 | tr ',' '\n' | tr '-' '\t' | \
awk '{print "chr10\t" $1 "\t" $2}' > a.bed
"""



################################################################################
##### 1. read DGs, find blockers, build models, compute correlation
def SEblockers(asmdict,fadict,DGdiv,phastCons,samples,
               indices,flank,loopmax,DGext,topDGs,stdmin):
    """
    one function to get DGs, build basepair models, and extract conservation.
    add new info dict to asmdict step by step.
    asmdict from readrMATS() in struclib.py
    input: asmdict[AScoord]['ASdict'] = {cell:[incm,incs,sum_IJC_SJC],...}
    output update:
    asmdict[AScoord]['covdicti'] = covdicti, intron cov dict for normalization
    asmdict[AScoord]['duplexes3X'] = duplexes3X
    asmdict[AScoord]['duplexes5'] = duplexes5
    covdicti = {sample:count} #here count is the sum of HNRNPC asso and indep
    duplexes3X/5 = [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,
                   seq,struc,mfe,pairsg,label1,label2,cons1,cons2}

    loopmax=50, max len of loops that block core motifs, instead of two arms
    flank=40, extending from SS3 to include PPT and BP, i.e. SS3X region. 
    DGext=10, extending each arm on each side by 10nts or until two arms merge.
    topDGs=10 or int(1E9), top DGs by cov, the most likely blockers.
    mininc=0; maxinc=1, parameters used in readrMATS(), not here. 
    stdmin=0, inc level std min, 0.1 is reasonable based on prior tests
    
    """
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    def divcomma(string): return [int(i) for i in string.split(',')]
    DGdivf = pyBigWig.open(DGdiv,'r'); phastf = pyBigWig.open(phastCons,'r')
    for AScoord in asmdict: #process each AS event one by one. 
        chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn = AScoord
        if chrom not in DGdivf.chroms(): continue
        ASdict = asmdict[AScoord]['ASdict']
        std = np.std([i[0] for i in ASdict.values()])
        if std <= stdmin: continue #take most dramatic ones for initial tests. 





        #if chrom != 'chr10' or ES>92479200 or EE<92478000: continue #SE in IDE
        if chrom != 'chr21': continue #test run time on one chrom. 
        #if chrom != 'chr19' or ES>1652615 or EE<1609292: continue #TCF3





        
        #A1. get gap1 cov in the two introns as norm standards.
        i1 = DGdivf.entries(chrom,EEup,ES); i1 = [] if not i1 else i1
        i2 = DGdivf.entries(chrom,EE,ESdn); i2 = [] if not i2 else i2
        entriesi = set(i1+i2) #merging left and right side introns
        if not entriesi: continue
        covdicti = {s:0 for s in samples} #cov dict for two introns combined
        for entry in entriesi:
            s1,e1,string = entry; data = string.split()
            sizes,starts = divcomma(data[7]),divcomma(data[8])
            DGl,DGr = [s1,s1+sizes[0]],[e1-sizes[1],e1]
            if overlap(*DGl,EEup,ES) or overlap(*DGr,EEup,ES) or \
               overlap(*DGl,EE,ESdn) or overlap(*DGr,EE,ESdn):
                for s in samples: covdicti[s] += int(data[indices[s]])

        #A2. get gap1 cov in SE, take geomean w/ intronic cov as norm standard
        covdicte = {s:0 for s in samples} #cov dict for the SE only
        entries = DGdivf.entries(chrom,ES,EE);
        if not entries: continue
        for entry in entries:
            s1,e1,string = entry; data = string.split()
            sizes,starts = divcomma(data[7]),divcomma(data[8])
            DGl,DGr = [s1,s1+sizes[0]],[e1-sizes[1],e1]
            if overlap(*DGl,ES,EE) or overlap(*DGr,ES,EE):
                for s in samples: covdicte[s] += int(data[indices[s]])
        for s in samples: covdicti[s] = (covdicti[s]*covdicte[s])**0.5
        asmdict[AScoord]['covdicti'] = covdicti

        #B. get all DGs overlapping SS5, a single nt position
        SS5 = (chrom,EE-1,EE) if strand=='+' else (chrom,ES,ES+1)
        entries5 = DGdivf.entries(*SS5); DGs5 = []; #all DGs blocking SS5
        duplexes5 = [] #all duplexes blocking SS5, 8 elements per DG now: 
        #[DGID,l1,l2,r1,r2,sumc1,sumc2,covdict]
        if entries5:             
            for entry in entries5: #DGs
                s1,e1,string = entry; data = string.split(); DGID = data[0]
                sizes,starts = divcomma(data[7]),divcomma(data[8])                
                if s1<SS5[1]<=s1+sizes[0] or e1-sizes[1]<=SS5[1]<e1 or \
                   s1+sizes[0]<SS5[1]<e1-sizes[1] and \
                   e1-sizes[1]-(s1+sizes[0])<=loopmax: #overlaps SS5
                    covdict = {k:data[indices[k]] for k in indices}
                    c1 = [int(data[indices[k]]) for k in indices if asso in k]
                    c2 = [int(data[indices[k]]) for k in indices if indep in k]
                    nums = [s1,s1+sizes[0],e1-sizes[1],e1]
                    DGs5.append([DGID]+nums+[sum(c1),sum(c2),covdict])

            #C. get models for top DGs that overlap SS5
            if DGs5:
                #12 elements per DG now, c1/c2: HNRNPCasso/indep. 
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
                    if ssl1<SS5[1]<ssl2 or ssr1<SS5[1]<ssr2 or \
                       ssl2<=SS5[1]<=ssr1 and ssr1-ssl2<=loopmax:#overlap duplex
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

                        #it is not fair to compare intronic vs. exonic cons.
                        #label them in this order: SS5, I(ntron), E(xon).
                        labels = []
                        if overlap(ssl1,ssl2,*SS5[1:3]): labels.append('SS5')
                        elif overlap(ssl1,ssl2,EEup,ES) or \
                             overlap(ssl1,ssl2,EE,ESdn): labels.append('I')
                        else: labels.append('E')
                        if overlap(ssr1,ssr2,*SS5[1:3]): labels.append('SS5')
                        elif overlap(ssr1,ssr2,EEup,ES) or \
                             overlap(ssr1,ssr2,EE,ESdn): labels.append('I')
                        else: labels.append('E')
                        duplexes5[i] = duplex+labels+[cons1,cons2]
        asmdict[AScoord]['duplexes5'] = duplexes5
        

            
        #E. get all DGs overlapping BP+PPT+SS3, e.g. a 40nt region
        SS3X = (chrom,ES-flank,ES) if strand=='+' else (chrom,EE,EE+flank)
        entries3X = DGdivf.entries(*SS3X); DGs3X = []
        duplexes3X = [] #all duplexes blocking SS3X. 8 elements per DG now: 
        #[DGID,l1,l2,r1,r2,sumc1,sumc2,covdict]
        if entries3X:
            for entry in entries3X: #DGs
                s1,e1,string = entry; data = string.split(); DGID = data[0];
                sizes,starts = divcomma(data[7]),divcomma(data[8])
                if overlap(s1,s1+sizes[0],*SS3X[1:]) or \
                   overlap(e1-sizes[1],e1,*SS3X[1:]): #overlaps SS3X
                    covdict = {k:data[indices[k]] for k in indices}
                    c1 = [int(data[indices[k]]) for k in indices if asso in k]
                    c2 = [int(data[indices[k]]) for k in indices if indep in k]
                    nums = [s1,s1+sizes[0],e1-sizes[1],e1]
                    DGs3X.append([DGID]+nums+[sum(c1),sum(c2),covdict])
            
            #F. get models for top DGs overlapping SS3X
            if DGs3X:
                #from 8 to 12 elements now per DG, c1/c2: HNRNPCasso/indep
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
                    if overlap(ssl1,ssl2,*SS3X[1:]) or \
                       overlap(ssr1,ssr2,*SS3X[1:]):
                        duplexes3X.append(DG+list(ss)+[pairsg])

                #G. get conservation for top duplexes (not DG) that overlap SS3X
                if duplexes3X:
                    #16 elements now per DG: [DGID,l1,l2,r1,r2,sumc1,sumc2,
                    #covdict,seq,struc,mfe,pairsg,label1,label2,cons1,cons2]
                    for i in range(len(duplexes3X)):
                        duplex = duplexes3X[i]; pairsg = duplex[-1]
                        if pairsg[0] == [0,0]:
                            duplexes3X[i] = duplex+['NA','NA',0,0]; continue
                        ssl1,ssl2,ssr1,ssr2 = sorted(pairsg[0]+pairsg[-1])
                        cons1 = phastf.stats(chrom,ssl1,ssl2)[0]
                        cons2 = phastf.stats(chrom,ssr1,ssr2)[0]
                        armcons = 0 if not (cons1 and cons2) else cons1*cons2 #?
                        labels = []
                        if overlap(ssl1,ssl2,*SS3X[1:3]): labels.append('SS3X')
                        elif overlap(ssl1,ssl2,EEup,ES): labels.append('I')
                        elif overlap(ssl1,ssl2,EE,ESdn): labels.append('I')
                        else: labels.append('E')
                        if overlap(ssr1,ssr2,*SS3X[1:3]): labels.append('SS3X')
                        elif overlap(ssr1,ssr2,EEup,ES): labels.append('I')
                        elif overlap(ssr1,ssr2,EE,ESdn): labels.append('I')
                        else: labels.append('E')
                        duplexes3X[i] = duplex+labels+[cons1,cons2]
        asmdict[AScoord]['duplexes3X'] = duplexes3X
        
    DGdivf.close(); phastf.close(); return asmdict #worked well. 
################################################################################





################################################################################
def exportasm(asmdict,samples,outfile): #export to a txt file
    """
    asmdict[AScoord]:
    {'ASdict':ASdict,'covdicti':covdicti,'duplexes3X':[],'duplexes5':[]}
    ASdict: {cell:[incm,incs,sum_IJC_SJC],...}
    output 21 columns in 5 groups: DGout, gap1out, ssout, ASout, consout
    gap1cv1/2 may not provide sufficient info for sorting across AS events,
    due to the variable coverage among introns. 
    """
    f = open(outfile,'w'); nDGs = 0
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    header = ['DGID','DGcoord','asso,indep','gap1asso,gap1indep','gap1total',
              'gap1cv1','gap1cv2','gap1block','gap1introns','seq','struc',
              'mfe','basepairs','SE','blocksite','incstd','incdict',
              'phastCons1','phastCons2','phastConsmean1','phastConsmean2']
    f.write('\t'.join(header)+'\n')
    for AScoord in asmdict:
        if 'covdicti' not in asmdict[AScoord]: continue
        keys = 'ASdict','covdicti','duplexes3X','duplexes5'
        ASdict,covdicti,duplexes3X,duplexes5 =[asmdict[AScoord][k]for k in keys]
        incstd = str(np.std([ASdict[c][0] for c in ASdict]))
        AScoordstr = ','.join([str(i) for i in AScoord])
        #note: not every AS event has inclusion data from all cell lines
        ASdictstr = ';'.join(['{}={}'.format(c,ASdict[c]) for c in ASdict])
        covistr = ';'.join(['{}={}'.format(s,covdicti[s]) for s in covdicti])
        #covistr, cov dict intron string, 30 samples, all DGs in two introns
        exps = 'HNRNPCassomul,HNRNPCindepmul' #experiments

        duplexes = duplexes3X+duplexes5
        for i in range(len(duplexes)):
            nDGs+=1; blocksite = 'SS3X' if i<len(duplexes3X) else 'SS5'
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
    f.close(); return nDGs
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
        sys.exit("python spliceblockSE.py fasta DGdiv phastCons suffix")
    #4 types of input: fasta,DGdiv,phastCons,rmatsfolders
    #rmatsfolders = ["astro_vs_neuron_outBAM_igvgtf_rmats"] #quick test
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

    #D. read the splicing data from rMATS output. ~34s for 10 SE comparisons. 
    asmdict = readrMATS(cwd,rmatsfolders,'SE',mininc,maxinc);
    print(str(datetime.datetime.today())+'\tFinished making splicing dict')
    
    #E. read DGs, get blockers, build models, get phastCons, modify asmdict
    #new data are added to the asmdict. 
    data = asmdict,fadict,DGdiv,phastCons,samples,indices
    params = flank,loopmax,DGext,topDGs,stdmin
    asmdict = SEblockers(*data,*params); #print(asmdict)
    print(str(datetime.datetime.today())+'\tFinished processing all data')

    #F. export data, 21 columns. Examples are analyzed later.
    nDGs = exportasm(asmdict,samples,outtxt)
    print(str(datetime.datetime.today())+'\tFinished exporting blockers')
    #for k in asmdict: print(asmdict[k]['ASdict']); break; sys.exit()
    #key=('chrY','+',ES,EE,ESup,EEup,ESdn,EEdn)
    #'ASdict': {'astro':[1,0,585],'neuron':[0.998,0.0026,234]},...}
    
    print(); filters = (topDGs,mininc,maxinc,stdmin)
    print("Filters: topDGs={},mininc={},maxinc={},stdmin={}".format(*filters))
    print("Other: flank={},loopmax={},DGext={}".format(flank,loopmax,DGext))
    print("Number of AS events passing filter: {}".format(len(asmdict)))
    print('Exported DG blockers: {}\n'.format(nDGs))
    
    #G. computes correlation and plot correlation scatters 
    #see details in a separate script: splicevsDG.py
    #summarize data using this script: ... 
################################################################################







"""
##### hg38 all 5 samples. How much time? ~3-4 hours locally.
Next test cutoff at 0.95, do we need cutoff for top5 or more?
What if we set cutoff to 1? All potential blockers, many of which may not work.
Set topDGs=1E9 to remove this limit, how many are there in total. 
Include all the 7 cell lines
Filters: topDGs=5,maxinc=0.9,stdmin=0; Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 54770
Exported DG blockers: 453601, 10220 genes, 50954 SE events. 23510065 gap1 reads
f=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_sorted_Birch_T20_\
sorted_DGdiv_iPS5lines_SEblockers.txt #1.0G
We identified base pairing regions that may cover these regions. Highest ones,
dozens of blockers per core motif. Lowest DG coverage regions, 0 or 1 per motif. 

##### chr22 for 5 samples: 372 sec. 
Filters: topDGs=5,maxinc=0.9,stdmin=0; Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 54770; Exported DG blockers: 10546

ATE1 MXE example is ranked very high on the list, based on MFE. 
PLKHA1, good example, not the best gap1 coverage.
TCF7L2, very high coverage of gap1, strong AS, likely SE and MXE. high cons. 
crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_\
sorted_Birch_T20_sorted_DGdiv_SEblockers_chr22.txt
52M in size. 23910 DGs. 



##### Manual inspection of output from chr10, setting stdmin=0.2
PFKFB3, chr10:6,225,764-6,231,759, strong duplex, blocking branch point
not tested correlation across cell types yet. not very obvious
PFKP, chr10:3,119,314-3,131,305, strong blocking, not very clear in diff. cells

Analysis of TCF3 identified a few DGs with (near-)significant correlation with
inclusion levels at the switching region. 
crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_sorted_Birch_T20_\
sorted_DGdiv_iPS5lines_SEblockers_TCF3.txt
We need to use both statistical tests to find more examples. 


##### TCF3:
Filters: topDGs=5,mininc=0,maxinc=0.9,stdmin=0
Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 54770
Exported DG blockers: 109, 262k file size

2025-10-10 14:33:01.153071	Finished processing all data
2025-10-10 14:33:01.401641	Finished exporting blockers
Filters: topDGs=1000000000,mininc=0,maxinc=1,stdmin=0
Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 127364
Exported DG blockers: 6092, just for one gene, way too many, 13M file size. 
"""


