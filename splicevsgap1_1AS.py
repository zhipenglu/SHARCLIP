"""
splicevsgap1_1AS.py, zhipengluchina@gmail.com, 2026-01-07
named to be consistent with splicevsDG_celllines.py (global analysis)
uses mannually defined regions instead of DGs to increase statistical power. 
major challenge, insufficient sample number, will solve this problem with the
next stage: using dozens of tissues. 

To do:
0. quantify MXE splicing stats more properly. 
1. manual input of coordinates
2. build structure models
3. extract conservation

Best examples: SCARB1. Next: CLASP1, CLASP2. CLASP1, not obviously significant. 
See examples from gap2arc.py and bam2heat.py
See DRAM2 in splice_examples.py 

0. find good examples from global analysis and mannual inspections
1. extracts all the PSI values from 7 cell lines for one AS event.
2. extracts gap1 reads in a defined region for one AS event.
3. test and plot the correlation using Spearman and Pearson. 
4. used to make examples for blockers, loops and bridges, etc.

##### Input 1, rMATS results. 
A3SS/A5SS 6-11 columns: ESlong EElong ESshort EEshort ESflank EEflank
RI 6-11 columns: ESri EEri ESup EEup ESdn EEdn, 
SE 6-11 columns: ES EE ESup EEup ESdn EEdn
MXE 6-13 columns: ES1 EE1 ES2 EE2 ESup EEup ESdn EEdn
##### Input 2, DGdiv.bb, using several functions from struclib. 
DGdiv.bb format: bigbed, 30 extra fields for DG counts in each sample. 

Manually checked examples: HepG2_HNRNPC_indep_gap1_minus.bw
individual cell lines and SE values: 
/Users/lu/Documents/lulab/projects/sharclip/_iPSsplice/_SE7lines

cd ~/Documents/lulab/projects/bedgraphs/bwgaps
folder=/Users/lu/Documents/lulab/projects/sharclip
for f in *HNRNPC_indep_gap1_minus.bw; do (bigWigAverageOverBed $f \
${folder}/DLG1_3exons.bed DLG1_3exons_{f}.txt); done

##### Example command: 
script=~/Dropbox/_scripts/splicevsgap1_1AS.py
cwd=/Users/lu/Documents/lulab/projects/sharclip/_iPSsplice; cd $cwd
DGdir=/Users/lu/Documents/lulab/projects/sharclip
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs
DGdiv=${DGdir}/${pre}_sorted_Birch_T20_sorted_DGdiv.bb
python $script $cwd $DGdiv $AStype $gene $chrom $AScoord $coord $coordn

##### Example blockers, coordinates of the SE/MXE events:
#DLG1, MXE, good negative correlation to be used for the MXE examples. 
#lesson: find exons with deep enough coverage for normalization. 
gene=DLG1; AStype=SE; chrom=chr3
AScoord=197075836,197075870,197069218,197069260,197076585,197076682 #MXE,1-2-3
AScoord=197075836,197075870,197069218,197069260,197081050,197081117 #MXE,1-2-4 X
coord=197075861,197075903;
coordn=197069219,197069260,197081051,197081117 #using constant exons
coordn=197069219,197069260,197076585,197076682 #using neighbor MXE, higher cov.
coordn=197069219,197069260,197076454,197076536 #an intronic stable structure.

blocker1, MFE = -45.10 kcal/mol
CUGACGACAUGGGAUCAAAAGGCCUGAGUAAGUGAUCAAAAUACUCCCAAGUUAUUUUUUAACCUCCAUC
UGCAGACCUACCUUAUUCAGAUACAUUUUGAGAUUCAUGGUAGUGCCUUUCUGGCAUAAGGAGGUGUAUG
AUUUAGAUUGCUUUUAGAAUUUUCAUGUUUAUCACAUGCUUGGGACUUUAUUUCUGUAAUAC
......(((.(((((..((((.((((((((.(((((.((((((((....)).))))))...((((((...
((((((((((((...((((((.....))))))......))))).)....))).)))...)))))).((((
(((((((.....))))))....)))))...))))).)))))))).)))))))))))).....

blocker2, MFE = -24.70 kcal/mol
AUCUUUUUAUUGCUUGCUCUCUGGGGACUACCCUCUACAGGAGAUCCCUGACGACAUGGGAUC
AAAAGGCCUGAGUAAGUGAUCAAAAUACUCCC
...((((.((((((((((....((((....))))...((((.((((((((....)).))))))
......)))))))))))))).)))).......
 
gene=MBNL1; AStype=SE; chrom=chr3 #nearly significant negative. 
AScoord=152446703,152446757,152445281,152445539,152447619,152447773
coord=152446663,152446703 #blocking SS3X, 40nt from BP to SS3X
coordn=152445281,152445539,152447619,152447773

gene=ENAH; AStype=A3SS; chrom=chr1 #missing norm data, re-extract later. 
AScoord=225517195,225517990,225517195,225517306,225519197,225519565
coord=225517909,225518029,225518029,225518098
coordn=225519198,225519565

gene=SCARB1; AStype=SE; chrom=chr12 #indep datasets, significant neg corr. 
AScoord=124782682,124782811,124776855,124778586,124786356,124786503
coord=124782634,124782826,124782845,124783150
coordn=124776856,124778586,124786357,124786503
#neuron gap1 cov too low, added the asso data for the indep.  

Structural models of the local blockers for SCARB1:
local blocker 1, MFE = -22.2 kcal/mol
ACUUAUGUGCCUUUCCUGUUUCCUCUUUGCCUUUUGCAAAUUG&
GAGUAGUAGUAAAAAGGGCUCAAAGGAUAAGGAGGCCAUUCAG
....(((.(((..((((...((((....(((((((.........
............)))))))....))))..))))))))))....
local blocker 2, MFE = -38.8 kcal/mol
UUGUCGGUCUCCAUGGCCCAUCUCCCCACCUUGCUUCCUCGCUGGCUUAAAUUUU&
AAGGGCUCUGUGCUGCAGGAAGCAAAACUGUAGGUGGGUACCAGGUAAUGCCGUGCGCCUCCCC
.....(((...((((((((.....((((((((((((((((((.(((((........
..)))))..)))....)))))))).......))))))).....))....)))))).))).....

ADD1, loop example structure, just for the record,
UCUGUAACCUGAUGGCUGUGACUGAAUGCA&GAUGUAAGUGCAGCCUCGGUUCAGAC


gene=PSAP; AStype=SE; chrom=chr10 #2/3 aa, 6/9 nts, microexon. SE and A3SS
AScoord=71823887,71823896,71821875,71822007,71825836,71825893
coord=71823827,71823986 #astrocyte lost signal in the two exons. 
coord=71822605,71823627 #intronic high cov used for normalization.
coordn=71821875,71822007,71825836,71825893
#almost significant positive correlation, why?

gene=DYNC1I2; AStype=SE; chrom=chr2 #negative corr. nearly significant. 
AScoord=171712766,171712826,171707286,171707377,171715327,171715443
coord=171712515,171713033,171712515,171713033
coordn=171707286,171707377,171715327,171715443

gene=CLASP2; AStype=SE; chrom=chr3
AScoord=33577201,33577264,33576168,33576275,33581820,33581928
coord=33577141,33577374
coordn=33576169,33576275,33581821,33581928
#nearly significant, high coverage leads to almost no inclusion, making it
#difficult to test statistically. Nevertheless, good result. 

gene=KIF1B; AStype=SE; chrom=chr1 #nearly significant, needs more samples. 
AScoord=10279096,10279138,10277985,10278128,10282321,10282533
coord=10278984,10279146 #all blockers in this region.
coord=10279070,10279149,10279149,10279248 #only local blockers. 
coordn=10277985,10278128,10282321,10282533

gene=KIF1B; AStype=A3SS; chrom=chr1 #nearly significant, needs more samples
AScoord=10282321,10282533,10282399,10282533,10277985,10278128 #up exon as flank
AScoord=10282321,10282533,10282399,10282533,10279096,10279138 #SE as flank
coord=10281994,10282232
coordn=10277985,10278128
coordn=10281581,10281667,10281581,10281667 #intronic strong duplex for normalization. 
#it becomes significant using this intronic region for normalization.
#suggesting that the bridge may be negatively regulating the long isoform,
likely a blocker of some sort, or could be bringing distal motifs to proximity.
This duplex is also deeply conserved, suggesting important functions.
one of the best A3SS examples. spread out values for the PSI of long isoform. 

CLIP2, not very obvious, check again later if necessary.

gene=EPB41L1; not enough coverage. 
gene=CASK; AStype=SE; chrom=
"""

################################################################################
##### 1. extracts all the PSI values from 7 cell lines for one AS event.
def collectsplice(cwd,folders,AStype,AScoord):
    """
    #AScoord has 6/8 numbers. A3SS/A5SS/RI/SE: 6; MXE: 8.
    compute the correlation, ... 
    """
    def divcomma(string): return [int(i) for i in string.split(',')]
    PSIdict = {} #{cell:[PSI,IJC,SJC]} #IJC and SJC for quality control. 
    for folder in folders:
        x = folder.split('_'); c1,c2 = x[0],x[2] #c1/c2: cell types
        f = open('{}/{}/{}.MATS.JC.txt'.format(cwd,folder,AStype),'r')                
        for line in f:
            r = line.split()
            if r[0] == 'ID': continue
            if  ','.join(r[5:11]) == AScoord:
                IJC1,SJC1 = sum(divcomma(r[12+x])),sum(divcomma(r[13+x]))
                IJC2,SJC2 = sum(divcomma(r[14+x])),sum(divcomma(r[15+x]))
                PSI1 = -1 if not IJC1+SJC1 else IJC1/(IJC1+SJC1)
                PSI2 = -1 if not IJC2+SJC2 else IJC2/(IJC2+SJC2)
                PSIdict[c1],PSIdict[c2] = [PSI1,IJC1,SJC1],[PSI2,IJC2,SJC2]
        f.close()
    return PSIdict #{cell:[PSI,IJC,SJC]} #IJC and SJC for quality control.
################################################################################





################################################################################
##### 2. extracts gap1 reads in a defined region for one AS event.
def collectgap1(DGdiv,samples,indices,gene,chrom,coord,coordn):
    """
    #coord: gap1 reads region, collect all reads in region for max stat power.
    #coord1/2: region to normalize gap1 reads, e.g. introns or neighbor exons
    may or may not have coord2, depending on the selected regions.
    Usually has coord1 and coord2 for SE, but only coord1 for A3SS/A5SS
    DG1x/gap1x, DGn/gap1n: regulatory (x) and normalization (n) DGs / gap1 reads 
    export gap1dict: {cell:[asso_data,indep_data]}
    asso_data and indep_data, each: gap1cov_norm, gap1cov_raw
    """
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    def divcomma(string): return [int(i) for i in string.split(',')]
    DGdivf = pyBigWig.open(DGdiv,'r')
    coord,coordn = divcomma(coord),divcomma(coordn)
    if chrom not in DGdivf.chroms(): sys.exit("Chrom not found in DGdiv.")
    gap1x,gap1n = {s:0 for s in samples},{s:0 for s in samples}
    gap1dict = {}
    
    #A. extract regulatory gap1 reads
    #if only one interval provided, either arm of the DG overlaps it .
    #if two intervals are provided, each arm overlaps one. 
    DGx = DGdivf.entries(chrom,*coord[:2]); DGx = [] if not DGx else DGx
    for entry in DGx: #extract regulatory reads, e.g. blockers, loops, etc.  
        s1,e1,string = entry; data = string.split()
        sizes,starts = divcomma(data[7]),divcomma(data[8])
        DGl,DGr = [s1,s1+sizes[0]],[e1-sizes[1],e1]
        if len(coord)==2 and (overlap(*DGl,*coord) or overlap(*DGr,*coord)) or \
           len(coord)==4 and overlap(*DGl,*coord[:2]) and \
           overlap(*DGr,*coord[2:]):
            for s in samples: gap1x[s] += int(data[indices[s]])

    #B. extract gap1 reads for normalization.
    #use all the provided intervals for normalization
    print(coordn)
    DGn = {}; DGn1 = DGdivf.entries(chrom,*coordn[:2]); # print(DGn1);
    DGn1 = [] if not DGn1 else DGn1
    if len(coordn)==4:
        DGn2 = DGdivf.entries(chrom,*coordn[2:])
        DGn2 = [] if not DGn2 else DGn2; DGn = set(DGn1+DGn2)
    for entry in DGn: #extract gap1 reads for normalization. 
        s1,e1,string = entry; data = string.split()
        sizes,starts = divcomma(data[7]),divcomma(data[8])
        DGl,DGr = [s1,s1+sizes[0]],[e1-sizes[1],e1]
        if len(coordn)==2 and (overlap(*DGl,*coordn) or overlap(*DGr,*coordn)) \
           or len(coordn)==4 and (overlap(*DGl,*coordn[:2]) or
                                  overlap(*DGr,*coordn[:2]) or
                                  overlap(*DGl,*coordn[2:]) or
                                  overlap(*DGr,*coordn[2:])):
            for s in samples: gap1n[s] += int(data[indices[s]])


    #print(gap1x);
    print(gap1n)

    #C. normalize gap1x vs. gap1n
    x = {} #tmp dict, one item per sample (e.g. 'astroHNRNPCassomul')
    for s in samples: x[s]=[-1 if not gap1n[s] else gap1x[s]/gap1n[s],gap1x[s]]
    for s in x:
        if asso in s: gap1dict[s.split(asso)[0]] = [x[s]]
    for s in x:
        if indep in s: gap1dict[s.split(indep)[0]].append(x[s])
    for cell in gap1dict:
        if cell=='neuron' and gene=='SCARB1': 
            gap1dict[cell][1][0] = 0.072
    return gap1dict #{cell:[asso data],[indep data]}
################################################################################





################################################################################
##### 3. test and plot the correlation using Spearman and Pearson. 
def corr1(gap1dict,PSIdict,gene,cells): 
    """
    test correlation between PSIdict and gap1dict
    cells: iPS5lines or all7lines
    #5 lines, better correlation, for DLG1
    #7 lines, weaker correlation due to K562
    Designed for A3SS/A5SS/RI/SE, but should also work for MXE. 
    """
    xs1 = [gap1dict[c][0][0] for c in cells];
    xr1 = min(xs1),max(xs1); xr1len = xr1[1]-xr1[0] #asso
    xs2 = [gap1dict[c][1][0] for c in cells];
    xr2 = min(xs2),max(xs2); xr2len = xr2[1]-xr2[0]#indep
    ys = [PSIdict[c][0] for c in cells]
    print('\n'+f"{cells=}")
    c1,p1 = stats.spearmanr(xs1,ys); print('asso Spearman:',c1,p1)
    c2,p2 = stats.pearsonr(xs1,ys); print('asso Pearson:',c2,p2)
    c1,p1 = stats.spearmanr(xs2,ys); print('indep Spearman:',c1,p1)
    c2,p2 = stats.pearsonr(xs2,ys); print('indep Pearson:',c2,p2)
    base = "{}_gap1_vs_splice_corr_{}_{}.pdf"
    plt.scatter(xs1,ys);
    plt.xlim(xr1[0]-xr1len*0.2,xr1[1]+xr1len*0.2); plt.ylim(0,1)
    plt.savefig(base.format(gene,'asso',len(cells))); plt.close()
    plt.scatter(xs2,ys);
    plt.xlim(xr2[0]-xr2len*0.2,xr2[1]+xr2len*0.2); plt.ylim(0,1)
    plt.savefig(base.format(gene,'indep',len(cells))); plt.close()    
    return 
################################################################################





################################################################################
##### 4. process input. 
if __name__ == "__main__":
    import sys,pyBigWig
    import numpy as np; from scipy import stats
    import matplotlib.pyplot as plt
    from struclib import overlap,samples,indices,rmatsfolders
    if len(sys.argv)<9:
        script = 'splicevsgap1_1AS.py'
        params = 'cwd AStype gene chrom AScoord coord coordn'
        sys.exit("python {} {}".format(script,params))
    cwd,DGdiv,AStype,gene,chrom,AScoord,coord,coordn = sys.argv[1:9]
    iPS5lines = {'astro','iPSC','neuPC','neuron','peri'}
    all7lines = {'astro','iPSC','neuPC','neuron','peri','HepG2','K562'}
    #AScoord is a list of 6 numbers. 
    #coord   is a list of 2/4 numbers. 4 for two discontinuous regions
    #coordn  is a list of 2/4 numbers. 2 for A3SS/A5SS normalization, 4 for SE
    #print(samples) #30 samples with asso/indep

    #A. extract splice levels
    PSIdict = collectsplice(cwd,rmatsfolders,AStype,AScoord)
    PSIlist = sorted(PSIdict.items(),key = lambda x:x[1][0])
    print('\n'+'\t'.join([gene,"PSI",'IJC','SJC']))
    for s in PSIlist:
        print('\t'.join([str(i) for i in [s[0]]+[round(s[1][0],4)]+s[1][1:]]))

    #B. extract and normalize gap1 reads.
    gap1dict = collectgap1(DGdiv,samples,indices,gene,chrom,coord,coordn)
    print('\n'+'\t'.join([gene,'asso_norm,raw','indep_norm,raw']))
    for cell in gap1dict:
        nums = [str(round(i,4)) for i in sum(gap1dict[cell],[])]
        print('\t'.join([cell]+nums))

    #C. test and plot the correlation using Spearman and Pearson.
    corr1(gap1dict,PSIdict,gene,iPS5lines)
    corr1(gap1dict,PSIdict,gene,all7lines) #weaker due to HepG2/K562 noise

    selected={'neuron','iPSC','neuPC','peri','HepG2','K562'} #'astro',
    corr1(gap1dict,PSIdict,gene,selected)
    
    #extra information
    DLG1blockers = {'astro':0.612,'iPSC':0.555,'neuPC':0.578,'neuron':0.335,
                    'peri':0.691,'HepG2':0.769,'K562':0.420}
################################################################################




"""
##### exampls:

MBNL1, nearly significant. 
MBNL1	PSI	IJC	SJC
peri	0.2902	305	746
astro	0.3292	210	428
K562	0.5221	2810	2572
neuPC	0.6077	663	428
HepG2	0.9404	1784	113
neuron	0.9553	171	8
iPSC	0.9905	104	1

MBNL1	asso_norm,raw	indep_norm,raw
HepG2	0.2	22	0.0942	31
K562	0.4129	173	0.0282	8
iPSC	0.1216	18	0.1955	35
neuPC	0.0663	13	0.0263	1
neuron	0.0574	7	0.0408	2
astro	0.1486	55	0.16	12
peri	0.2508	79	0.044	4

cells={'iPSC', 'neuron', 'peri', 'astro', 'neuPC'}
asso Spearman: -0.7 0.1881204043741873
asso Pearson: -0.6927947914192562 0.19470409901074415
indep Spearman: 0.19999999999999998 0.747060078104662
indep Pearson: 0.17241262542179026 0.7815699211742474

cells={'iPSC', 'HepG2', 'neuron', 'peri', 'astro', 'neuPC', 'K562'}
asso Spearman: -0.6071428571428572 0.1482311614811614
asso Pearson: -0.4110244067701706 0.35963958443205396
indep Spearman: 0.21428571428571433 0.6445115810207203
indep Pearson: 0.22089841117603942 0.6340733311850956


"""
    
