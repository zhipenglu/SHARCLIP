"""
spliceblockRI.py, 2025-10-02, zhipengluchina@gmail.com
see parameters within script. Replacing previous two versions spliceblockSE.py.
Rewritten on 2025-10-04 from spliceblockSE.py. To be tested carefully. 
See Monteuuis 2019 for a review of RI: The changing paradigm of intron
retention: regulation, ramifications and recipes

Changes for spliceblockRI.py
1. only one intron considered ...
2. SS3X and SS5 coordinates are different from before. 

##### Basic procedure. 
1. read rMATS output AS data from 7 cell lines, readrMATS() 
2. read DG from DGdiv.bb
3. build bp models for top DGs that block core motifs
4. extract conservation from phastCons
5. compute correlation among inc levels, gap1 reads, and structure stability. 

##### Major function:
RIblockers(asmdict,fadict,DGdiv,phastCons,exp,indices)
A. Get DGs in two regions: SS5 and BP-PPT-SS3
B. To compare genes in one sample, use MFE, no need for normalization? TO DO.
C. To compare one AS event across samples, use cov when MFE is OK. 

#conditions for blocking SS5:
A. For DGs, overlaps either arm, or in the loop and loop size<=50nts
B. For duplex models: overlaps either arm, or in the loop and loop size<=50nts

#conditions for blocking SS3 extended region (SS3 + 40nts, SS3X):
A. For DGs, overlaps either arm
B. For duplex models: overlaps either arm. What about loop region? Modify later.

##### Find good examples of RI.
BICD2, strong RI, big diff across cell types. structure data so so.

FNBP4, chr11:47731300-47732700, 1.4kb, deep conservation, complex retention
patterns, a super strong structure blocking one 3' ss. both SE and RI. 
CUUCACCGAGGUUAAUCAGGAGUAGUUGCUUCCAGAGGUCCUGUAAAGCGCACAGAGCU\
CUCGCAUGCGCUUAACCCCCAAGCCCGCCCUACUCCGUGCAACACUCUGGAGAGUAA
(((((..((((((...(((((((((..((((...(.((....((.(((((((..((...\
.))...))))))).))))).)))).....))))))).)).))).)))))))).....

PRPF4B, chr6:4,060,317-4,062,229, strong ~3kb conservation of 3'UTR retention,
similar to the FNBP4, strong structure ~200nts range. Full length: 1289nts. 
PRPF4B. chr6 + 4060412 4061700 4060412 4061147 4061654 4061700
GUACACUGAAGGUAAGCUAAACCAUCAACAUCUCUGGUGUUUUAAGA&
CAUUUAAAAUUGGAAGAGAACGCGCUUGAUGGAUAGAGCGCCUUCAGUGUACUGU
(((((((((((((..(((...((((((((.(((((...(((((((...
...)))))))....)))))....).)))))))....))))))))))))))))...

HLTF, chr3:149,031,180-149,031,734, strong duplex, SS3 in the loop region.
ENDOV, significant diff among cells. not conserved. chr17:80429568-80430499
NFKB2, clear diff among cells, chr10:102,395,669-102,396,042, gap1 data so so

##### Example command
fasta=~/igv/genomes/hg38/hg38.fa
fasta=~/Documents/lulab/projects/hg38/fasta/chr11.fa
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs
DGdiv=${pre}_sorted_Birch_T20_sorted_DGdiv.bb
phastCons=~/Documents/lulab/projects/hg38/conservation/hg38.phastCons100way.bw
suffix=_all7lines_RIblockers_all.txt
python ~/Dropbox/_scripts/spliceblockRI.py $fasta $DGdiv $phastCons $suffix &


"""



################################################################################
##### 1. read DGs, find blockers, build models, compute correlation
def RIblockers(asmdict,fadict,DGdiv,phastCons,samples,
               indices,flank,loopmax,DGext,topDGs,stdmin):
    """
    one function to get DGs, build models, and extract conservation.
    add new info dict to asmdict step by step.
    asmdict from readrMATS(), from struclib.py 
    input: asmdict[AScoord]['ASdict'] = {cell:[incm,incs,sum_IJC_SJC],...}
    output update:
    asmdict[AScoord]['covdicti'] = covdicti
    asmdict[AScoord]['duplexes3X'] = duplexes3X
    asmdict[AScoord]['duplexes5'] = duplexes5
    covdicti = {sample:count} #here count is the sum of HNRNPC asso and indep
    duplexes3X/5 = [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,
                   seq,struc,mfe,pairsg,label1,label2,cons1,cons2}

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
        chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn = AScoord #ES==ESup; EE==EEdn.
        if chrom not in DGdivf.chroms(): continue
        ASdict = asmdict[AScoord]['ASdict']
        std = np.std([i[0] for i in ASdict.values()])
        if std <= stdmin: continue #take most dramatic ones for initial tests.

        #PRPF4B. chr6 + 4060412 4061700 4060412 4061147 4061654 4061700
        #if chrom != 'chr11' or EE<47731300 or ES>47732700: continue #FNBP4






        #A. get gap1 cov in the retained intron as norm standards.
        #Should we use the retained or surrounding introns for normalization?        
        entriesi = DGdivf.entries(chrom,EEup,ESdn) #retained intron
        if not entriesi: continue
        covdicti = {s:0 for s in samples} #cov dict for two introns combined
        for entry in entriesi:
            s1,e1,string = entry; data = string.split()
            sizes,starts = divcomma(data[7]),divcomma(data[8])
            DGl,DGr = [s1,s1+sizes[0]],[e1-sizes[1],e1]
            if overlap(*DGl,EEup,ESdn) or overlap(*DGr,EEup,ESdn):
                for s in samples: covdicti[s] += int(data[indices[s]])
        asmdict[AScoord]['covdicti'] = covdicti
        
        #B. get all DGs overlapping SS5, a single nt position
        #for RI, only changed the SS5 coords the first line below. 
        SS5 = (chrom,EEup-1,EEup) if strand=='+' else (chrom,ESdn,ESdn+1)
        entries5 = DGdivf.entries(*SS5); DGs5 = []; #all DGs blocking SS5
        duplexes5 = [] #all duplexes blocking SS5.
        #8 elements per DG now: [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict]
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
            #did not change anything from SE for this section
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
                    if ssl1<SS5[1]<ssl2 or ssr1<SS5[1]<ssr2 or \
                       ssl2<=SS5[1]<=ssr1 and ssr1-ssl2<=loopmax:
                        duplexes5.append(DG+list(ss)+[pairsg])
                
                #D. get conservation for top duplexes (not DGs) overlapping SS5
                if duplexes5:
                    #16 elements per DG now, [DGID,l1,l2,r1,r2,sumc1,sumc2,
                    #covdict,seq,struc,mfe,pairsg]
                    for i in range(len(duplexes5)):
                        duplex = duplexes5[i]; pairsg = duplex[-1]
                        if pairsg[0] == [0,0]:
                            duplexes5[i] = duplex+['NA','NA',0,0]; continue
                        ssl1,ssl2,ssr1,ssr2 = sorted(pairsg[0]+pairsg[-1])                        
                        cons1 = phastf.stats(chrom,ssl1,ssl2)[0]
                        cons2 = phastf.stats(chrom,ssr1,ssr2)[0]

                        #it is not fair to compare intronic vs. exonic cons.
                        #label them properly: I(ntron), E(xon), SS5.
                        labels = []
                        if overlap(ssl1,ssl2,*SS5[1:3]): labels.append('SS5')
                        elif overlap(ssl1,ssl2,EEup,ESdn): labels.append('I')
                        else: labels.append('E')
                        if overlap(ssr1,ssr2,*SS5[1:3]): labels.append('SS5')
                        elif overlap(ssr1,ssr2,EEup,ESdn): labels.append('I')
                        else: labels.append('E')
                        duplexes5[i] = duplex+labels+[cons1,cons2]
        asmdict[AScoord]['duplexes5'] = duplexes5
        
        #E. get all DGs overlapping BP+PPT+SS3, e.g. a 40nt region
        SS3X = (chrom,ESdn-flank,ESdn)if strand=='+' else(chrom,EEup,EEup+flank)
        entries3X = DGdivf.entries(*SS3X); DGs3X = []
        duplexes3X = [] #all DGs blocking SS3X with some structure.
        #8 elements per DGnow: [DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,seq]
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
                #12 elements per DG now, c1/c2: HNRNPCasso/indep
                #[DGID,l1,l2,r1,r2,sumc1,sumc2,covdict,struc,mfe,pairsg]
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
                    #16 elements per DG now, [DGID,l1,l2,r1,r2,sumc1,sumc2,
                    #covdict,seq,struc,mfe,pairsg,label1,label2,cons1,cons2]
                    for i in range(len(duplexes3X)):
                        duplex = duplexes3X[i]; pairsg = duplex[-1]
                        if pairsg[0] == [0,0]:
                            duplexes3X[i] = duplex+['NA','NA',0,0]; continue
                        ssl1,ssl2,ssr1,ssr2 = sorted(pairsg[0]+pairsg[-1])
                        cons1 = phastf.stats(chrom,ssl1,ssl2,exact=True)[0]
                        cons2 = phastf.stats(chrom,ssr1,ssr2,exact=True)[0]
                        labels = []                        
                        if overlap(ssl1,ssl2,*SS3X[1:3]): labels.append('SS3X')
                        elif overlap(ssl1,ssl2,EEup,ESdn): labels.append('I')
                        else: labels.append('E')
                        if overlap(ssr1,ssr2,*SS3X[1:3]): labels.append('SS3X')
                        elif overlap(ssr1,ssr2,EEup,ESdn): labels.append('I')
                        else: labels.append('E')
                        duplexes3X[i] = duplex+labels+[cons1,cons2]
        asmdict[AScoord]['duplexes3X'] = duplexes3X
    DGdivf.close(); phastf.close(); return asmdict #worked well. 
################################################################################





################################################################################
def exportasm(asmdict,samples,outfile): #export to a txt file
    """
    asmdict[AScoord]: {'ASdict':ASdict,'duplexes3X':[],'duplexes5':[]}
    ASdict: {cell:[incm,incs,sum_IJC_SJC],...}
    output 21 columns in 5 groups: DGout, gap1out, ssout, ASout, consout
    gap1std1/2 do not provide sufficient info for sorting across AS events,
    due to the variable coverage among introns. 
    """
    f = open(outfile,'w'); nDGs = 0
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    header = ['DGID','DGcoord','asso,indep','gap1asso,gap1indep','gap1total',
              'gap1cv1','gap1cv2','gap1block','gap1introns','seq','struc',
              'mfe','basepairs','RI','blocksite','incstd','incdict',
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
            nDGs+=1
            blocksite = 'SS3X' if i<len(duplexes3X) else 'SS5'
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
    import numpy as np; from matplotlib import pyplot as plt
    from struclib import readrMATS,readfa,ssmodel,db2pairs,overlap,makesamples
    from struclib import cwd,rmatsfolders,samples,indices

    if len(sys.argv)<4:sys.exit("python spliceblockRI.py fasta DGdiv phastCons")
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

    #D. read the splicing data from rMATS output. ~34s for 10 SE comparisons. 
    #5 cell lines, n=186102 for SE. After filtering with maxinc=0.98, n=40005.
    asmdict = readrMATS(cwd,rmatsfolders,'RI',mininc,maxinc);
    print(str(datetime.datetime.today())+'\tFinished making splicing dict')
    
    #E. read DGs, find blockers, build models, extract phastCons
    #new data are added to the asmdict. 
    data = asmdict,fadict,DGdiv,phastCons,samples,indices
    params = flank,loopmax,DGext,topDGs,stdmin
    asmdict = RIblockers(*data,*params); #print(asmdict)
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
    #see details in a separate script: spliceblockcorr.py
################################################################################





"""
Tested on example FNBP4, strongly conserved RI with variable inclusion across
cell lines. Very strong structure. 

2025-10-14 09:07:41.726935	Starting analysis
2025-10-14 09:07:58.177111	Finished reading fasta
2025-10-14 09:07:59.745161	Finished making splicing dict
2025-10-14 09:38:19.035427	Finished processing all data
2025-10-14 09:38:22.382563	Finished exporting blockers
Filters: topDGs=10,mininc=0,maxinc=1,stdmin=0
Other: flank=40,loopmax=50,DGext=10
Number of AS events passing filter: 6042
Exported DG blockers: 73578

"""
