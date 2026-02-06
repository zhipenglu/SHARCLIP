"""
DG2model.py; zhipengluchina@gmail.com, 2025-10-10
To be used with gap2arc.py and bam2heat.py to visualize data.
Can be used to visually identify strong loops, blockers, bridges and switches
gap2arc.py and bam2heat.py quantify connection strength better than DG2model.py 
Application to transcriptome is not necessary at the moment
Input data: fasta and DGdiv.bb.
Note: automatically built models from DGs may be quite different from manually
built ones. DGs even if extended by 10nts on each end may still not be long
enough to cover the important regions. TCF7L2 is a good example where the
manually build models are much longer, ~100bps.

##### Which samples to use?
1. HNRNPCindepmul is better for most purposes.
2. Other RBPs work well for regions enriched in other RBPs. 
3. Primarily uses asso and indep as groups, instead of individual cell types.

##### Focusing on the potential switch for PKM. Similar filter for other genes.
awk '$1=="track" || $2>72200900 && $2<72201010 && $3<72202340' \
PKM_bps_indep_covlim0.01.bed > PKM_bps_indep_covlim0.01_switch2.bed

##### Example command:
cd ~/Documents/lulab/projects/sharclip; script=~/Dropbox/_scripts/DG2model.py;
mamba activate crssant; f=~/igv/genomes/hg38/hg38.fa
pre=crssant_all_sharclip_eclip_sharc_paris_all_gap1_IGV_mRNAs_sorted
DGdiv=${pre}_Birch_T20_sorted_DGdiv.bb
d=~/Documents/lulab/projects/hg38/fasta/

f=${d}chr3.fa; c='HNRNPCindep FOXP1 chr3 70970000 9000 0.1'
f=${d}chr9.fa; c='HNRNPCindep DNM1 chr9 128222000 13000 0.1'
f=${d}chr1.fa; c='HNRNPCindep PTBP2 chr1 96804500 3000 0.05' #UCE/SE.
f=${d}chrY.fa; c='HNRNPCindep UTY chrY 13357800 2200 0.1'
f=${d}chr9.fa; c='HNRNPCindep TPM2 chr9 35684400 1000 0.2'
f=${d}chr11.fa; c='HNRNPCindep DPF2 chr11 65343700 2700 0.1'
f=${d}chr12.fa; c='HNRNPCindep KIF21A chr12 39326000 5000 0.02' #also RBFOX2
f=${d}chr15.fa; c='HNRNPCindep PKM chr15 72200000 7000 0.02'
f=${d}chr10.fa; c='HNRNPCindep IDE chr10 92478200 2200 0.1'
f=${d}chr19.fa; c='HNRNPCindep TCF3 chr19 1611500 4500 0.02'
f=${d}chr10.fa; c='HNRNPCindep TCF7L2 chr10 113158000 10000 0.05'

python $script $f $DGdiv $c


"""


################################################################################
##### 1. read DG information, get coverage for a set of samples
def readDGdiv(DGdiv,exp,indices,chrom,start,end):
    #exp: a full name or part of it, e.g. astroHNRNPCindepmul or indep
    #indices: {sample_name:column_id}
    #output DGdict: {DGID:[chrom,strand,l1,l2,r1,r2,cov]}
    def divcomma(string): return [int(i) for i in string.split(',')]
    f = pyBigWig.open(DGdiv,'r'); genome = f.chroms()
    entries = f.entries(chrom,start,end); DGdict = {}
    print(chrom,start,end)
    if not entries: return DGdict
    for entry in entries: #each entry is one DG in bed12+30 format.
        s1,e1,string = entry; data = string.split()
        if s1 < start or e1 > end: continue
        DGID,strand = data[0],data[2]#print([k for k in indices if exp in k])
        sizes,starts = divcomma(data[7]),divcomma(data[8])
        cov = sum([int(data[indices[k]]) for k in indices if exp in k])
        if cov: DGdict[DGID] = [chrom,strand,s1,s1+sizes[0],e1-sizes[1],e1,cov]
    return DGdict
################################################################################





################################################################################
def makemodels(fadict,DGs,outfile):
    #make bp models and export to bed file. 
    bedout = []
    for DG in DGs[::-1]:
        chrom,strand,l1,l2,r1,r2,cov = DG
        seq,struc,mfe = ssmodel(chrom,[l1,l2,r1,r2],strand,fadict)
        pairs = db2pairs(struc); sep = list(struc).index('&')
        for pair in sorted(pairs):
            pair[0] = pair[0]+l1 if pair[0]<sep else pair[0]+r1
            pair[1] = pair[1]+l1 if pair[1]<sep else pair[1]+r1
            c = str(int((1-cov/DGs[0][6])*255)); c = ','.join([c,c,c])
            extra = '.\t1\t.\t0\t0\t'+c #not needed here. 
            bedout.append([chrom,str(pair[0]),str(pair[1])])
    string = '\n'.join(['\t'.join(i) for i in bedout])+'\n'
    f = open(outfile,'w'); f.write('track graphType=arc\n'+string); f.close()
    return
################################################################################





################################################################################
##### 5. process files. 
if __name__ == "__main__":
    import os; os.environ['OPENBLAS_NUM_THREADS'] = '1' #override a bug on CARC
    import sys, pyBigWig; import numpy as np; from datetime import datetime
    from matplotlib import pyplot as plt
    from struclib import readfa,ssmodel,db2pairs,makesamples

    if len(sys.argv) < 7:
        sys.exit("DG2model.py fasta DGdiv prefix sample chrom start size lim")
    fasta,DGdiv,exp,gene,chrom,start,size,covlim = sys.argv[1:9]
    start,size,covlim = int(start),int(size),float(covlim); end = start+size
    outfile = 'bps_{}_{}_{}_{}_{}_covlim{}.bed'.format(sys.argv[3:9])
    samples = makesamples() #sample name, e.g. HepG2HNRNPCassomul
    indices = {samples[i]:i+9 for i in range(len(samples))} #e.g. 
    
    #A. get genome fasta file, #faidx samtools?
    fadict = readfa(fasta) 

    #B. read DG info from DGdiv.bb
    DGdict = readDGdiv(DGdiv,exp,indices,chrom,start,end)
    covs = sorted([i[6] for i in DGdict.values()])
    #plt.plot(covs); plt.savefig("a.pdf") #check distribution

    #C. make models
    DGs1 = sorted(list(DGdict.values()),key=lambda x:x[6],reverse=1);
    DGs2 = [DG for DG in DGs1 if DG[6]>=covlim*DGs1[0][6]];
    makemodels(fadict,DGs2,outfile)
    print("Number of DGs before and after filtering", len(DGs1),len(DGs2))
################################################################################





"""
#1492 DGs in the TCF3 region, max cov = 341,
#DGs with cov >=10% max: 94; with cov >=2% max: 408
    

"""
