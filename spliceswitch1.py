"""
spliceswitch1.py, zhipengluchina@gmail.com, 2025-09-11
replacing the older version spliceswitch.py,
minor changes: using exp to control analysis of HNRNPCasso vs. HNRNPCindep data
include structure modeling and switch discovery. 

From MXE and SHARCLIP data, identify structural switches that control them.
See spliceswitchscan.py for simple alternatives to find alt. blockers and loops,
which can be done on individual MXE clusters or all of them together.
See script: DG2model.py for structure modeling.
Here we also build helices and alt. conformations. 
ATE1 was not included in Hatje 2017 list, likely due to a few DI reads.
~/Desktop/sharclip/_iPSsplice/astro_vs_neuron_outBAM_igvgtf_rmats/
MXE.MATS.JC.txt, 38013 records, much more than the Hatje 2017 annotation.
The other annotation by Martinez Gomez 2021 reported ~230 MXE clusters, which
are much more reliable.

##### note on 2025-10-16, reorganize data, export to a file, similar to the
spliceblock*.py output format. Likely very different due to two DGs per switch. 

##### Lift Hatje 2017 annotations from hg19 to hg38. See also reorg_Hatje.py. 
awk -F '\t' '{print $3 "\t" $5 "\t" $6 "\t" $9 "," $10 "\t" int($14) "\t" $7}' \
Hatje2017_s2.txt > Hatje2017_s2.bed
chain=~/Documents/lulab/projects/hg38/hg19ToHg38.over.chain 
liftOver Hatje2017_s2.bed $chain Hatje2017_s2_hg38.bed Hatje2017_s2_unmapped

##### Input:
1. MXE data from rMATS, how to store them? Start with examples, ATE1/TCF3. 
2. MXE annotation from Hatje 2017.
3. associated biology and disease. Hatje2017_s7_SNP.xls, ~30 genes, likely more
4. DG data from DGdiv.bb, focus on 5 iPSC-derived cell lines or all 7 cell lines
5. fasta, for structure prediction, try next. 
6. phastCons and phyloP, additional layer of ranking.

##### Basic procedure: 
1. read MXE output from rMATS, readrMATS()
2. intersect with Hatje 2017 annotations, filterMXE()
4. get all DGs from DGdiv, rank by coverage, readDGs()
5. predict structures for top ranked ones, gethelices()
6. identify alt. conformations, findalt() #not done yet. 

##### Expected output, lists, summary figures, and examples:
0. for ~1400 MXEs, we have quantified inc levels for a small subset, ~560
1. the Martinez Gomez 2021 list of 236 clusters is much more reliable.
2. a list of DGs blocking or looping the alt. exons, with cov across samples
3. structure models, esp. stable helices, defined by MFE
4. correlation of DG cov with inclusion levels.
5. intersect with conserved regions, 100way datasets
6. ~100 MXE clusters with good evidence: DGs, structures, conservation
7. rank by alt. helix strength, strong differences in DG cov, and conservation
7. to be plotted: structure stability of primary blocking and looping DGs
8. nominate new therapeutic targets for genes where MXE is linked to diseases. 

##### select a subset of DGs? What criteria?
1. arms must overlap introns, at least partially. 
1. must span two introns OR overlap splicing motifs
3. determine the mechanism (block or loop) after identifying functional ones. 

##### What are the reasonable analysis?
Visual inspections and statistical summary suggest <100 confidence MXE clusters.
Each cluster associate with hundreds of DGs even after filtering
for ones that are close to core splicing motifs. The dominant ones may not be
all that obvious. We need to decide on proper criteria for calling switches.

##### Output file formats, see example from spliceblockSE.py
Normalization can be performed between alt. conformations, so no need to extract
data for entire introns. How do we identify pairs and report them?
Label them as pairs in column 22? For example, altDGs=1 for all DGs in 1 group.
Some DGs may be part of multiple groups so multiple altDGs values are needed. 

##### Example command:
cd ~/Documents/lulab/projects/sharclip/; mamba activate crssant
script=~/Dropbox/_scripts/spliceswitch1.py
MXEanno=~/Dropbox/_2025_SHARCLIP/Hatje2017/Hatje2017_s2_hg38.bed
MXEall=~/Dropbox/_2025_SHARCLIP/Hatje2017/Hatje2017_s2.txt
fasta=~/Documents/lulab/projects/hg38/fasta/chr19.fa
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_sorted
DGdiv=${pre}_Birch_T20_sorted_DGdiv.bb
phylop=~/Documents/lulab/projects/hg38/conservation/hg38.phyloP100way.bw
phastCons=~/Documents/lulab/projects/hg38/conservation/hg38.phastCons100way.bw

time python $script $MXEanno $MXEall $fasta $DGdiv $phastCons

"""


################################################################################
##### 1. compute reading frames, i.e. is exon length a multiple of 3?
"""
assuming all the MXEs in each cluster have the same in-frame status.
could be wrong for non-homologous ones. enable multiple options.
Do we need to check the conservation down to flies? not really
"""
def readframe(infile):
    #Hatje2017_s2.txt
    f = open(infile,'r'); framedict = {}
    for line in f:
        if line.split()[0] == 'Accession': continue
        record = line.split('\t'); frame,name,ID = record[7:10];
        frame = [int(i) for i in frame.split('-')] #'x-y', with x,y in [0,1,2]
        framedict[(name,ID)] = 1 if frame[0]-frame[1]==0 else 0 #1:in-frame
    f.close()
    framedict[('ACTN4','0')] = 0; framedict[('ATE1','0')] = 1
    return framedict
################################################################################





################################################################################
##### 2. Other information: U2U12 compatibility and disease annotations.
#disease1: MXE biased, from Hatje2017, Hatje2017_s7_SNP.xls
#disease2: other direct disease connections, manual search'
#Later we mannually checked all the diseases for the therapeutic development. 
U2U12 = {"MAPK8","MAPK9","MAPK10","MAPK14","CEP170","CRTC1"} #Hatje2017 Fig.S18
disease1 = {"TPM2":"Cap myopathy 2",
            "TTN":"myopathy 1G",
            "SNAP25":"Myasthenic",
            "MAPT":"PSP and FTD",
            "SLC25A3":"MPCD",
            "CACNA1C":"Long QT syndrome",
            "ADH1B":"Alcohol dependence",
            "FGFR3":"colon carcinoma",
            "SCN5A":"dilated cardiomyopathy",
            "FHL1":"EDMD6",
            "CACNA1C":"PFVF",
            "CACNA1D":"PASNA",
            "DNM2":"Charcot-Marie-Tooth",
            "TPM1":"PFHC",
            "FHL1":"Scapuloperoneal myopathy",
            "FGFR1":"Osteoglophonic dysplasia",
            "FGFR3":"Carcinoma of colon",
            "IQCB1":"SLSN5",
            "KCNQ2":"BFNS1",
            "EPCAM":"Lynch syndrome",
            "SNRNP200":"Retinitis pigmentosa",
            "SCN9A":"Primary erythromelalgia",
            "ACTN4":"FSGS1",
            "PNPLA6":"Boucher Neuhauser",
            "SCN8A":"epileptic encephalopathy",
            "MYL2":"FHC",
            "FAR1":"FAR1 disorder",
            "BSCL2":"Charcot-Marie-Tooth",
            "LDB3":"Dilated cardiomyopathy 1C",
            "FGFR2":"Crouzon syndrome",
            "GNAI3":"Auriculocondylar 1",
            "TCF3":"Burkitt's and immuno def",
            "PKM":"cancer"}
disease2 = {"SCN2A":"epileptic encephalopathies",
            "ACTN2":"myopathy",
            "P4HA1":"connective tissue disorder",
            "ANKRD36C":"iTTP",
            "ADAM23":"DEE",
            "DTNA":"muscular dystrophy",
            "COL6A2":"muscular dystrophy",
            "CACNA1A":"neuropathologies",
            "OGDH":"OGDH deficiency",
            "PDLIM3":"cardiomyopathy",
            "ACOX1":"ACOX1 deficiency",
            "NEB":"nemaline myopathy",
            "NGLY":"NGLY deficiency",
            "MTHFSD":"neurodev disorders",
            "OBSCN":"cardiomyopathy",
            "PDLIM7":"FMVP",
            "CDK11B":"syndactyly",
            "GTF2I":"Williams syndrome",
            "LDB3":"cardiomyopathy",
            "P4HA2":"connective tissue disease",
            "CACNA1E":"DEE",
            "AIFM1":"encephalomyopathy",
            "COL5A1":"classical Ehlers-Danlos",
            "SRGAP3":"intellectual disability",
            "CACNB1":"muscular disorder",
            "USH1C":"Usher 1C",
            "IFT43":"cranioectodermal dysplasia",
            "H2AFY":"Liebenberg",
            "DIS3":"multiple myoloma",
            "TPM3":"myopathy",
            "SLC39A14":"hypermanganesemia",
            "DLG3":"intellectual disability",
            "GRIA1":"epilepsy",
            "GRIA2":"NEDLIB",
            "GRIA3":"intellectual disability",
            "CAMK2D":"neuropath, cardiomyopathy",
            "UBE2E2":"diabetes",
            "TCTN1":"Joubert syndrome",
            "U2AF1":"myelodysplastic syndromes",
            "KCNMA1":"neuromuscular disorders",
            "MEF2C":"neurodev disorder, haploinsuf",
            "ITGA7":"muscular dystrophy",
            "TNNT3":"distal arthrogryposis",
            "TRDN":"TKOS, CVPT",
            "SCN3A":"neurological disorders",
            "DNM1":"encephalopathy",
            "NGLY1":"NGLY1-CDDG",
            "PHKB":"glycogen storage GSD IXb"}
################################################################################





################################################################################
##### 3. intersect asmdict with MXEanno to find MXEs with RNA-seq measurements
"""
asmdict from earlier function readrMATS(). MXEanno, input file from Hatje 2017
here the filters are Hatje2017 annotations, just to be comprehensive.
the Martinez Gomez 2021 annotations are much more reliable. In fact, our
filtering in the end essentially converged with the Gomez study.
the highest IJC+SJC criterion was used to simplify cases where the MXEs also
have other AS types, e.g. A3SS, A5SS, etc. which are quite common, e.g. TCF3.
"""
def filterMXE(asmdict,MXEanno,phastCons):
    #MXEs and clusters by Hatje 2017: 1392 628
    #MXEs and clusters after filtering: 571 453
    Hatje2017 = {} #step A. all annotated MXEs from Hatje 2017
    MXEdict = {}   #step B. dict of MXE data after intersecting with Hatje 2017
    C1 = {}        #step C. clusters of MXE using SE combo with highest IJC+SJC

    #A. Add missed MXE groups manually.
    f1,f2 = open(MXEanno,'r'),pyBigWig.open(phastCons,'r')
    a,b = ('chr19','+',38727945,38728026),('chr19','+',38728315,38728381)
    Hatje2017[a] = ('ACTN4','0',f2.stats(a[0],*a[2:4],exact=True)[0])
    Hatje2017[b] = ('ACTN4','0',f2.stats(b[0],*b[2:4],exact=True)[0])
    a,b = ('chr10','-',121898840,121898969),('chr10','-',121899865,121899994)
    Hatje2017[a] = ('ATE1','0',f2.stats(a[0],*a[2:4],exact=True)[0])
    Hatje2017[b] = ('ATE1','0',f2.stats(b[0],*b[2:4],exact=True)[0])

    #B. get Hatje 2017 annotation: chrom start end name,cluster score strand 
    #most MXEs are in the SE output, so mechanisms will still be covered
    for line in f1:
        record = line.split()
        if record[0] == "Accession": continue
        chrom,start,end,x,score,strand = record; name,ID = x.split(',')
        start,end,cluster = int(start),int(end),int(ID)
        cons = f2.stats(chrom,start,end,exact=True)[0]
        cons = 0 if not cons else cons
        Hatje2017[(chrom,strand,start-1,end)] = (name,ID,cons)
    nHatje2017 = len(set([i[1] for i in Hatje2017.values()]))
    f1.close(); f2.close()

    #C. intersect our MXE data from rMATS with Hatje2017 to make MXEdict
    #asmdict[SE]['ASdict'] = [{sample:[incm,incs,count],...},...] 
    #MXEdict[(name,cluster,cons,chrom,strand,ES,EE)] = [[SE,asmdict[SE]],...]
    #each MXEdict element corresponds to >=1 SE combinations from rMATS output
    #for example, TCF3 contains A3SS events in the MXE cluster. 
    for SE in asmdict:
        chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn = SE
        if SE[:4] not in Hatje2017: continue
        k = Hatje2017[SE[:4]]+SE[:4] #(name,cluster,cons,chrom,strand,ES,EE)
        MXEdict[k] = MXEdict.get(k,[]) + [[SE,asmdict[SE]['ASdict']]]

    #C. pick major exon combo because each Hatje MXE may have several SE combos
    #This does not deal with clusters with multiple MXEs e.g. ANKRD36C
    #most multiple-MXE clusters are likely artifacts or special types, and thus
    #can be ignored here. 
    #Combine MXEs to clusters: C1 {(name,cluster):[[SE,asmdict[SE],cons],...]}
    #Although ~300 clusters remain after filtering, most are likely not real.
    #many can be filtered out based on phastCons, but not necessary now. 
    for k in MXEdict: #options: choose SE with highest IJC+SJC counts
        x = sorted(MXEdict[k],key=lambda x:sum([i[-1]for i in x[1].values()]))
        C1[k[:2]] = C1.get(k[:2],[]) + [x[-1]+[k[2]]] #k[2] cons of MXE
    print('MXEs and clusters by Hatje 2017:', len(Hatje2017),nHatje2017)
    print("MXEs and clusters after filtering:", len(MXEdict),len(C1))
    return C1 #aka clusters1 {(name,cluster):[[SE,asmdict[SE],cons],...]}
################################################################################







################################################################################
##### 4. read DG data and only keep ones most likely regulating the MXEs
#how do we deal with multi-MXE clusters? Treat them separately.
#here the filters are len(exons)==4, and at least 1 DG per MXE cluster
def readDGs(DGdiv,clusters1,samples,indices,flank):
    """
    MXE clusters1: {(name,cluster):[[SE,asmdict[SE],cons],...]}
    asmdict[SE]['ASdict'] = [{sample:[incm,incs,count],...},...]
    SE = (chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn)
    multiple SEs are combined to define the MXE cluster.
    Low cov may not mean low confidence, e.g. for IDE, the lower side did not
    block the respective exon. It was biased towards one exon for all 7 lines.
    """
    def divcomma(string): return [int(i) for i in string.split(',')]
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul'
    f = pyBigWig.open(DGdiv,'r'); clusters2 = {} #keep clusters passing tests
    blockdict,loopdict = deepcopy(clusters1),deepcopy(clusters1)
    bridgedict = deepcopy(clusters1) #no need to associate with asmdict
    for key in clusters1: #key = (name,cluster)


        #X. this TCF3 cluster has 3 MXEs, two of which are A3SS. 
        if key==('TCF3','1424'):
            clusters1[key] = [clusters1[key][0]]+[clusters1[key][2]]



        #A. merge overlapping exons in each cluster of MXEs.
        #especially important for MXEs that include A3SS, A5SS, etc. 
        SEs = [i[0] for i in clusters1[key]];
        chrom,strand = clusters1[key][0][0][:2]
        intvls = sum([[list(i[2:4]),list(i[4:6]),list(i[6:8])] for i in SEs],[])
        exons = merge(intvls); start,end = exons[0][0],exons[-1][1]

        #B. get DG info for each cluster of MXEs. 
        #n=4 exons for 2-exon clusters1, E1,E2,E3,E4.
        #[1609291,1611849],[1612206,1612430],[1615284,1615520],[1615685,1615821]
        #E2e,E3e: extended on the 5' end by 40nts to cover BP and PPT
        if len(exons)!=4 and 'TCF3' not in key: continue
        if len(exons)<4: continue
        E1,E2,E3,E4 = exons; E2e,E3e = [E2[0]-flank,E2[1]],[E3[0]-flank,E3[1]]
        if strand == '-': E2e,E3e = [E2[0],E2[1]+flank],[E3[0],E3[1]+flank]
        
        entries = f.entries(chrom,start,end)
        if not entries: continue
        for entry in entries: #each entry is one DG
            s1,e1,string = entry; data = string.split();
            if s1<start or e1>end: continue
            sizes,starts = divcomma(data[7]),divcomma(data[8]);
            cov1,cov2 = 0,0; #total cov for asso/indep samples respectively 
            covdict = {} #cov for individual samples
            for s in samples:
                covdict[s] = int(data[indices[s]])
                if asso in s: cov1 += int(data[indices[s]])
                elif indep in s: cov2 += int(data[indices[s]])
            intvl = [chrom,strand,s1,s1+sizes[0],e1-sizes[1],e1] #DG intvl.
            DGinfo = intvl+[cov1,cov2,covdict]
            
            #C. get all DGs that potentially regulate MXEs. 
            #1. two midpoints span an extended exon: exon+40nts (to BP), t1
            #2. at least one midpoint must be intronic, t2 and t3
            ml,mr = s1+sizes[0]/2,e1-sizes[1]/2 #left/right arm midpoints.
            t1 = [ml<E2e[0]<mr,ml<E2e[1]<mr,ml<E3e[0]<mr,ml<E3e[1]<mr]#cross
            t2 = [E1[1]<ml<E2e[0],E2e[1]<ml<E3e[0],E3e[1]<ml<E4[0]]
            t3 = [E1[1]<mr<E2e[0],E2e[1]<mr<E3e[0],E3e[1]<mr<E4[0]] 
            if any(t1) and (any(t2) or any(t3)):
                if key not in clusters2: clusters2[key]=[clusters1[key],[]]
                clusters2[key][1].append(DGinfo)

            #D. get blockers for two MXEs, t1 for MXE1, t2 for MXE2
            #worked well for correlation analysis
            ex = 0 #high value of exclusion zone makes it worse. Stick to 0.
            arm1,arm2 = (s1,s1+sizes[0]),(e1-sizes[1],e1)
            t1a = overlap(*E2e,*arm1) and E2e[1]+ex<e1-sizes[1]<e1<E3e[0]-ex
            t1b = overlap(*E2e,*arm2) and E1[1]+ex<s1<s1+sizes[0]<E2e[0]-ex
            t2a = overlap(*E3e,*arm1) and E3e[1]+ex<e1-sizes[1]<e1<E4[0]-ex
            t2b = overlap(*E3e,*arm2) and E2e[1]+ex<s1<s1+sizes[0]<E3e[0]-ex
            if t1a or t1b: blockdict[key][0].append(DGinfo)
            elif t2a or t2b: blockdict[key][1].append(DGinfo) 
            #E. get only potential loops. l1/l2 for looping left/right MXEs.
            #worked well for correlation analysis
            l1 = E1[1]<s1<s1+sizes[0]<E2e[0]<E2e[1]<e1-sizes[1]<e1<E3e[0] 
            l2 = E2e[1]<s1<s1+sizes[0]<E3e[0]<E3e[1]<e1-sizes[1]<e1<E4[0]
            if l1: loopdict[key][0].append(DGinfo)
            elif l2: loopdict[key][1].append(DGinfo)
            #F. get potential bridges.
            if E2e[1]<s1 and e1<E3e[0]: bridgedict[key][0].append(DGinfo)
    f.close(); return clusters2, blockdict, loopdict, bridgedict
################################################################################





################################################################################
##### 3. dist distribution between MXEs. Short dist suggests steric hindrance
def MXEbridge(clusters2):
    #this function is only used for informational purposes. 
    for key in clusters2:
        name,ID = key; chrom,strand = clusters2[key][0][0][:2]
        SEs = [i[0] for i in clusters2[key][0]]
        intvls = sum([[list(i[2:4]),list(i[4:6]),list(i[6:8])] for i in SEs],[])
        exons = merge(intvls); start,end = exons[0][0],exons[-1][1]
        if len(exons)<4: continue
        #for E2,E3 in combinations(exons[1:-1],2):
        #    print('\t'.join([name,str(ID),str(E3[0]-E2[1])]))
    return 
################################################################################





################################################################################
##### 4. Compute all correlation and export all mechanisms here:
"""
Perform this before the structure modeling step, separate blockers and loops. 
MXE clusters1: {(name,cluster):[[SE,asmdict[SE],cons,covdict],...]}
SE = (chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn)
asmdict[SE]['ASdict'] = [{sample:[incm,incs,count],...},...]
covdict[sample] = int

blockdict: {(name,cluster):[[SE,asmdict[SE],cons,DGinfo],...]}
DGinfo = [chrom,strand,l1,l2,r1,r2,cov1,cov2,covdict]

Exported information:
A.
B.
C.

D. blocker fraction of each DG
E. blocker intronic phastCons
G. loop fraction of each DG
H. loop intronic phastCons
I. linear and bridged distances between MXEs
J. dominant bridge phastCons
K. fraction of each bridge DG.
O. published ASOs taking advantage of the switching, ...
"""

def corrMXE(blockdict,loopdict,bridgedict,
            U2U12,framedict,disease1,disease2,phastCons,samples,exp):
    asso,indep = 'HNRNPCassomul','HNRNPCindepmul';
    f = pyBigWig.open(phastCons,'r'); corrdict = {}
    #blockdict, key=(name,ID), values as follows:
    #'U2U12'=0/1, spliceosome compatibility, no/yes
    #'frame'=0/1, reading frame compabilility, no/yes
    #"disease" = str, disease name
    #"MXEcons" = float, phaseCons of MXEs
    #"mutbias" = ..., bias of mutations between MXEs, not used now. 
    #"blockercorr" = [], blocker correlation between gap1 and inc levels
    #"blockercons" = [], mean phastCons of intronic arm, list or dict, DGs
    #"blockerfrac" = [], fraction of reads in each DG.
    #"loopcorr" = []
    #"loopcons" = []
    #"loopfrac" = []
    #"bridgedist" = [], dist between SS5 and BP before/after bridging
    #"bridgecons" = [], 
    #"bridgefrac" = []
    #additional information, e.g. structure stability? top DG in each mech.
    
    for k in blockdict: #key = (gene_name,clusterID)
        if ('TCF3','1424')==k: blockdict[k]=[blockdict[k][0]]+[blockdict[k][2]] 
        if len(blockdict[k])!=2: continue
        corrdict[k] = {}
        SE1,asmdict1,SE2,asmdict2 = blockdict[k][0][:2]+blockdict[k][1][:2]
        chrom = SE1[0]; ES1,EE1,ES2,EE2 = SE1[2:4]+SE2[2:4] #MXE start/end
        #SE: (chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn)
        blockcov1,loopcov1 = blockdict[k][0][3:],loopdict[k][0][3:]
        blockcov2,loopcov2 = blockdict[k][1][3:],loopdict[k][1][3:]
        b1,b2,l1,l2 = blockcov1,blockcov2,loopcov1,loopcov2 #shorten names
        a = ('chr1','+',0,0,0,0,0,0) #filler info for empty DGs.
        b1 = [[*a,{'':0}]]if not b1 else b1; b2 = [[*a,{'':0}]]if not b2 else b2
        l1 = [[*a,{'':0}]]if not l1 else l1; l2 = [[*a,{'':0}]]if not l2 else l2
        #b1,b2,l1,l2: each is a list of dicts, where key = sample

        #A. get U2U12 information:
        corrdict[k]["U2U12"] = 1 if k[0] in U2U12 else 0
        #B. get reading frame info:
        corrdict[k]["frame"] = framedict[k]
        #C. get disease information:
        if k[0] in disease1: corrdict[k]["disease"] = (disease1[k[0]],"MXE")
        elif k[0] in disease2: corrdict[k]["disease"] = (disease2[k[0]],"other")
        else: corrdict[k]["disease"] = ('','')
        #D. get conservation of MXEs. 
        cons1 = f.stats(chrom,ES1,EE1,exact=True)[0]
        cons2 = f.stats(chrom,ES2,EE2,exact=True)[0]
        cons1 = 0 if not cons1 else cons1; cons2 = 0 if not cons2 else cons2
        corrdict[k]["MXEcons"] = [cons1,cons2]
        
        #E. sum all DGs or take the top DG for blocker and loop analysis.
        ss = samples; #print("b1[0]", k, b1[0],l1[0])
        maxb1={s:max(d[8].get(s,0)for d in b1)for s in ss if exp in s}
        maxb2={s:max(d[8].get(s,0)for d in b2)for s in ss if exp in s}
        sumb1={s:sum(d[8].get(s,0)for d in b1)for s in ss if exp in s}
        sumb2={s:sum(d[8].get(s,0)for d in b2)for s in ss if exp in s}
        maxl1={s:max(d[8].get(s,0)for d in l1)for s in ss if exp in s}
        maxl2={s:max(d[8].get(s,0)for d in l2)for s in ss if exp in s}
        suml1={s:sum(d[8].get(s,0)for d in l1)for s in ss if exp in s}
        suml2={s:sum(d[8].get(s,0)for d in l2)for s in ss if exp in s}
        
        #8 possible ratios: max vs. sum, blocker vs. loop, asso vs. indep.
        corrdict[k]["blockercorr"],corrdict[k]["loopcorr"] = [],[]
        for s in samples:
            #only use high cov junctions
            if exp not in s or s.split(exp)[0] not in asmdict1 or \
               s.split(exp)[0] not in asmdict2: continue
            incdata1 = asmdict1[s.split(exp)[0]]
            incdata2 = asmdict2[s.split(exp)[0]]
            if incdata1[2]+incdata2[2]<=0:
                corrdict[k]["blockercorr"].append(0)
                corrdict[k]["loopcorr"].append(0)
                continue #IJC+SJC counts
            
            incdiff = incdata1[0]-incdata2[0]
            #E1. only use high cov DGs, input data: max blocker asso or indep
            sub,add = maxb1[s]-maxb2[s],maxb1[s]+maxb2[s]
            if maxb1[s]+maxb2[s]<10: corrdict[k]["blockercorr"].append(0)
            else: corrdict[k]["blockercorr"].append(abs(incdiff-sub/add))
            #print('_'.join(k),incdiff,0 if not add else sub/add)
            #print to export for correlation computation in MXEscatter.py

            #E2. only use high cov DGs, input data: max loop asso or indep            
            sub,add = maxl1[s]-maxl2[s],maxl1[s]+maxl2[s]
            if maxl1[s]+maxl2[s]<10: corrdict[k]["loopcorr"].append(0)
            else: corrdict[k]["loopcorr"].append(abs(incdiff-sub/add))
            print('_'.join(k),incdiff,0 if not add else sub/add)
            #print to export for correlation computation in MXEscatter.py
            
            """
            #E3. use all DGs, input data: sum blocker asso or indep
            sub,add = sumb1[s]-sumb2[s],sumb1[s]+sumb2[s]
            if sumb1[s]+sumb2[s]<10: continue #gap1 reads 
            corrdict[k].append([incdiff,sub/add])
            #E4. use all DGs, input data: sum loop asso or indep            
            sub,add = suml1[s]-suml2[s],suml1[s]+suml2[s]
            if suml1[s]+suml2[s]<10: continue #gap1 reads 
            corrdict[k].append([incdiff,sub/add])
            """

        col = 6 if exp=="HNRNPCasso" else 7 #HNRNPCindep. Two possibilities
        #picking the asso or indep DGs for ranking purposes. 

        #F. blocker intronic arm conservation, using top 5 DGs, ranked based on
        #cov2, i.e., all the HNRNPCasso or HNRNPCindep samples. 
        #use the lower value to represent the intronic one, which makes sense.
        #it is extremely rare for exonic cons. to be lower than intronic cons.
        corrdict[k]['blockercons'] = []
        DGs = blockdict[k][0][3:]+blockdict[k][1][3:] #2 sets of DGs for 2 MXEs
        #print("blockdict[k][0][3]",blockdict[k][0][3]) #example DG
        DGs = sorted(DGs,key=lambda x:x[col]) #top 5 DGs based on cov1 or cov2
        for DGinfo in DGs[-5:]: #DGinfo = [l1,l2,r1,r2,cov1,cov2,covdict]
            chrom,strand,l1,l2,r1,r2,cov1,cov2,covdict = DGinfo
            cons1 = f.stats(chrom,l1,l2,exact=True)[0]
            cons2 = f.stats(chrom,r1,r2,exact=True)[0]
            cons1 = 0 if not cons1 else cons1; cons2 = 0 if not cons2 else cons2
            corrdict[k]['blockercons'].append(min(cons1,cons2))
        #G. blocker DG fraction of reads
        total = sum([i[col] for i in DGs])#using asso or indep only
        corrdict[k]['blockerfrac'] = [i[col]/total if total else 0 for i in DGs]


        #GG. find alt. blockers, i.e. blocker switches, 
        """
        Algorithm:
        for each set of blocker DGs linked to a group of MXEs (DGs list):
        there are very likely to be alt. blockers given homology between MXEs.
        In the IDE example, we find clear structural switches between two MXEs.
        They are all biased towards one side, clear evidence. What can we plot?
        We already have the abs(incdiff - gap1diff), a proxy for switching
        Next, identify the
        Individual examples are also shown. Identify the values clearly? 
        

        Check the gap2arc.py output files, 
        
        """




        
        #H. loop intronic arms conservation
        corrdict[k]['loopcons'] = []
        DGs = loopdict[k][0][3:]+loopdict[k][1][3:] #two sets of DGs. 
        #print("loopdict[k][0][3]",loopdict[k][0][3]) #example DG
        DGs = sorted(DGs,key=lambda x:x[col])[-5:] #top 5 DGs based on cov1/cov2
        for DGinfo in DGs: #DGinfo = [l1,l2,r1,r2,cov1,cov2,covdict]
            chrom,strand,l1,l2,r1,r2,cov1,cov2,covdict = DGinfo
            cons1 = f.stats(chrom,l1,l2,exact=True)[0]
            cons2 = f.stats(chrom,r1,r2,exact=True)[0]
            cons1 = 0 if not cons1 else cons1; cons2 = 0 if not cons2 else cons2
            corrdict[k]['loopcons'].append(np.mean([cons1,cons2]))
        #I. loop DG fraction of reads
        total = sum([i[col] for i in DGs])#using asso or indep only
        corrdict[k]['loopfrac'] = [i[col]/total if total else 0 for i in DGs]

        #II. find alt. loops, i.e. loop switches,
        """
        
        """



        

        #J. bridge arm conservation
        corrdict[k]['bridgecons'] = []; DGs = bridgedict[k][0][3:]
        if not DGs:
            corrdict[k]['bridgecons'] = [0]
            corrdict[k]['bridgefrac'] = [0]
            corrdict[k]['bridgedist'] = [ES2-EE1,ES2-EE1]; continue
        DGs = sorted(DGs,key=lambda x:x[col])[-5:] #top 5 DGs based on cov2 indep
        for DGinfo in DGs: #DGinfo = [l1,l2,r1,r2,cov1,cov2,covdict]
            chrom,strand,l1,l2,r1,r2,cov1,cov2,covdict = DGinfo
            cons1 = f.stats(chrom,l1,l2,exact=True)[0]
            cons2 = f.stats(chrom,r1,r2,exact=True)[0]
            cons1 = 0 if not cons1 else cons1; cons2 = 0 if not cons2 else cons2
            corrdict[k]['bridgecons'].append(np.mean([cons1,cons2]))
        #K. bridge fraction of reads for DG.
        total = sum([i[col] for i in DGs])#using asso or indep only
        corrdict[k]['bridgefrac'] = [i[col]/total if total else 0 for i in DGs]
        #L. inter-MXE distance before and after bridging: SS5-BP or SS5-SS3
        #which DG should we designate as the bridge? Good example COL25A1
        #most abundant one may be local and irrelavant
        DGs = [DG+[(DG[5]-DG[2])**2*DG[col]] for DG in DGs]
        DGtop = sorted(DGs,key=lambda x:x[-1])[-5:]
        bridge = max([i[5]-i[2] for i in DGtop])
        corrdict[k]["bridgedist"] = [abs(ES2-EE1),abs(ES2-EE1)-bridge]
        
    #plot all the mechanisms and related information.
    keys = ['U2U12','frame','disease','MXEcons','blockercorr','blockercons',
            'blockerfrac','loopcorr','loopcons','loopfrac','bridgedist',
            'bridgecons','bridgefrac']
    for k in corrdict:
        ks = '_'.join(k)
        U2U12,frame = str(corrdict[k]["U2U12"]),str(corrdict[k]["frame"])
        disease = '_'.join(corrdict[k]["disease"])
        MXEcons = '_'.join([str(i) for i in corrdict[k]['MXEcons']])
        blockercorr = '_'.join([str(i) for i in corrdict[k]['blockercorr']])
        blockercons = '_'.join([str(i) for i in corrdict[k]['blockercons']])
        blockerfrac = '_'.join([str(i) for i in corrdict[k]['blockerfrac']])
        loopcorr = '_'.join([str(i) for i in corrdict[k]['loopcorr']])
        loopcons = '_'.join([str(i) for i in corrdict[k]['loopcons']])
        loopfrac = '_'.join([str(i) for i in corrdict[k]['loopfrac']])
        bridgedist = '_'.join([str(i) for i in [max(corrdict[k]['bridgedist'])]])
        bridgecons = '_'.join([str(i) for i in corrdict[k]['bridgecons']])
        bridgefrac = '_'.join([str(i) for i in corrdict[k]['bridgefrac']])
        x = [ks,U2U12,frame,disease,MXEcons,blockercorr,blockercons,blockerfrac,
             loopcorr,loopcons,loopfrac,bridgedist,bridgecons,bridgefrac]
        #print('\t'.join(x))
    return 
################################################################################





################################################################################
##### 5. process input files. Major considerations in rewriting: 
"""
This process is similar to spliceblock*.py but functions are separated here.
Output differs from other mechanisms, e.g. blockers, bridges and switches. 
Correlation between inclusion levels and switches has not been tested yet.
We will optimize this script to produce better output files. Correlation
analysis will identify a small set. 
1. read all SE data from 7 cell lines, readrMATS()
2. filter using Hatje 2017 reference, filterMXE()
3. read DGs, keep potential regulators. Read diff sets for diff purposes.
4. build bp models, focus on a small set of genes, gethelices()
5. find alt. conformations, findalt()
6. TO DO: get conservation information. May not be necessary. 
summary of our SE/MXE data after intersection with Hatje 2017 annotations
Total SEs in 7 lines: 241387
MXEs and clusters by Hatje 2017: 1388 628
MXEs and clusters after filtering: 567 451
"""
if __name__ == "__main__":
    import sys,pyBigWig; from intervaltree import IntervalTree as IT
    from copy import deepcopy
    import numpy as np; from itertools import product,combinations
    from datetime import datetime; from matplotlib import pyplot as plt
    from struclib import readfa,ssmodel,db2pairs,merge,itermerge,altcheck
    from struclib import overlapmax,makesamples,readrMATS,overlap
    from struclib import cwd,samples,indices,rmatsfolders,iPS5lines,all7lines
    
    if len(sys.argv) < 7:
        params = "MXEanno MXEall fasta DGdiv phastCons exp"
        sys.exit("python spliceswitch1.py {}".format(params))
    MXEanno,MXEall,fasta,DGdiv,phastCons,exp = sys.argv[1:7]
    #exp = 'HNRNPCasso' or 'HNRNPCindep' 
    #set to asso/indep to produce output separately. 

    mininc,maxinc = 0,1 #not used here. 
    fadict = readfa(fasta)
    flank = 40 #50? extending the 3' splice sites to include BP and PPT. 
    covmin = 0.02 #relative to max of the region? 
    #rmatsfolders = ["astro_vs_neuron_outBAM_igvgtf_rmats"]
    #rmatsfolders = ["astro_vs_neuron_example"]

    #0. get reading frame information
    framedict = readframe(MXEall); #print(framedict); sys.exit()
    
    #1. read MXE data from rMATS output, using SE files, output asmdict. 
    #output asmdict: key=(chrom,strand,ES,EE,ESup,EEup,ESdn,EEdn) #coord.
    #values = {sample:[inclevel,std,sum_IJC_SJC]}.
    print(str(datetime.today())+'\tStarting analysis')
    asmdict = readrMATS(cwd,rmatsfolders,'SE',mininc,maxinc);
    print(str(datetime.today())+'\tFinished reading rMATS data')
    print("Total SEs in 7 lines:", len(asmdict)) #7 lines:241387. 

    #2. read MXE annotations from Hatje 2017, and filter MXE data. 
    #clusters1, {(name,cluster):[[SE,asmdict[SE],cons],...]}
    #test plotting of inclusion levels. Higher variations are more useful. 
    clusters1 = filterMXE(asmdict,MXEanno,phastCons); #plot_MXEstats(clusters1)
    print(str(datetime.today())+'\tFinished filtering MXE') #check example: 
    #for k in clusters1: print(len(clusters1[k]))

    #Double-checked script up to this step, as of 2025-10-22.
    #Currently focusing on spliceswitchscan.py as an alternative approach.
    #This script is useful for global analysis. 
    #3. read DGdiv.bb, filter to select likely DGs, update the clusters1 dict
    params = (DGdiv,clusters1,samples,indices,flank)
    clusters2,blockdict,loopdict,bridgedict = readDGs(*params)
    print(str(datetime.today())+'\tFinished reading DGdiv')
    print("MXE clusters after reading DGs:", len(clusters2))
    print("Filters here: 2-MXE clusters with at least one DG")

    #4. optional, check MXE distance. 
    #MXEbridge(clusters2)

    """
    for k in blockdict:
        if 'TCF3' in k: print(blockdict[k][0][3],loopdict[k][0][3]); break
    """
    
    #4. export correlation between blockers/loops and inclusion of MXEs. 
    #plot the correlation using a separate script scatter.py
    corrMXE(blockdict,loopdict,bridgedict,
            U2U12,framedict,disease1,disease2,phastCons,samples,exp)

    #plot the correlation output file MXEcorr.txt using script plotMXEcorr.py

    sys.exit()
    #4. build bp models, assemble helices, ignore protein-mediated ones now
    bedfile = "TCF3_bps_1611500_4500.bed"
    clusters2 = gethelices(fadict,clusters2,bedfile,covmin)
    print(str(datetime.today())+'\tFinished building bp models')
    #print("MXE clusters:", len(clusters2))
    #for c in clusters2: print(clusters2[c])
    #5. find alt conformations, export structure models and plot the summary
    findalt(clusters2,samples)
    print(str(datetime.today())+'\tFinished finding alt conformations')
################################################################################




"""
Annotations of some likely false positives from Hatje 2017:
NRG1, good MXEs, but missed, why? also as 3' end. among best alt loops.
MXE n=3: 5 listed but not necessarily real. 
CACNB2, likely 3 MXEs, but very different in sizes
PDLIM5, likely annotation errors
SLIT3, annotation errors
USH1C, only 2 after lifting to hg38. 
"""
        
