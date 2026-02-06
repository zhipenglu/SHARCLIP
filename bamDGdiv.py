"""
bamDGdiv.py, zhipengluchina@gmail.com, 2025-08-11
divide bam DGs to individual samples with counts
output format, bed12 plus 30 extra for bigbed storage and random access later. 
This output will make downstream analysis faster.
we can analyze correlation and plot them as well.
Input: all the SHARCLIP samples included:
5 RBPs, 7 cell lines, asso vs. indep, and drug treatment (not splicing data).

##### 30 fields after the first bed12: 
field 13: uint asso_count1;	"DG_count"
field 14: unit indep_count1;	"DG_count"
...
field 41: uint asso_count15;	"DG_count"
field 42: unit indep_count15;	"DG_count"

##### Basic procedure
1. convert bedpe of DG info to bed12, and then bigbed for random access
2. extract all reads for each gene from bam, count reads for each sample
3. extract the DG information from DG bigbed and export to a bed file.
4. convert the new bed file to bigbed using the SQL file for random access
5. use bigbed with divided DG counts for association with splicing

##### Prepare the bigbed file
script=~/Documents/scripts/duplex/bed12bedpe/bedpetobed12.py
genome=~/Documents/lulab/projects/hg38/hg38.chrom456.sizes
f=sharclip_ric_cric_eclip_karr_sharc_paris_all_gap1_mRNAs_Birch_T30_dg.bedpe
bedpe=${f/.bedpe/_GAPDH.bedpe} #for GAPDH
bed=sharclip_ric_cric_eclip_karr_sharc_paris_\
all_gap1_mRNAs_Birch_T30_dg_GAPDH.bed
python $script ${f/.bedpe/_GAPDH.bedpe} $bed
export LC_COLLATE=C
sort -k1,1 -k2,2n $bed |awk '$5=1000' >${bed/GAPDH/sorted_GAPDH}
bedToBigBed ... $genome ...
For bigbed format, the score must be between 0 and 1000: use 1000
The output bigbed file is ~3.3 to 4.2 G.

##### Example command for local test, output *_DGdiv.bed
cd /Users/lu/Desktop/sharclip; script=~/Dropbox/_scripts/bamDGdiv.py
genesfolder=~/Documents/lulab/projects/hg38/refGene/
genes=${genesfolder}hg38_refGene.mergedintvls_chr25.bed
genes=GAPDH.bed
bampre=sharclip_ric_cric_eclip_karr_sharc_paris_all_gap1_mRNAs_Birch_T30
bam=${bampre}_sorted_GAPDH.bam; DGbb=${bampre}_dg_sorted_GAPDH.bb
python $script $genes $bam $DGbb

##### command on CARC, bamDGdiv.py, for T20, T23 and T29:
script=/project/zhipengl_72/labshare/scripts/bamDGdiv.py
dir=/project/jianhuib_673/sharclip/crssant_all
pre=crssant_all_sharclip_eclip_sharc_paris
cd $dir/$pre/${pre}_all_gap1_IGV_mRNAs_crssant_self_structure
genes=/project/zhipengl_72/labshare/hg38/hg38_refGene.mergedintvls_chr25.bed
DGbb=${pre}_all_gap1_IGV_mRNAs_sorted_Birch_T29_dg.bb
bam=${pre}_all_gap1_IGV_mRNAs_sorted_Birch_T29_sorted.bam
#bam=GAPDH.bam; genes=GAPDH.bed
salloc --time=10:00:00 --mem=10G python $script $genes $bam $DGbb

##### command on CARC, sort and convert *_DGdiv.bed to *_DGdiv.bb:
bed=${pre}_all_gap1_IGV_mRNAs_sorted_Birch_T23_sorted_DGdiv.bed
#sorting took 16min, using script sortbed_DGdiv1.slurm. 
genome=/project/zhipengl_72/labshare/hg38/hg38.chrom456.sizes
bedToBigBed=/project2/zhipengl_72/labshare/kent/bedToBigBed; sql=DGdiv.as
salloc --time=10:00:00 --mem=10G $bedToBigBed -as=$sql -type=bed12+30 ${bed}3 \
$genome ${bed}4 &

##### Make sure all the chromosomes are included in the DGdiv.bed output !!!!!
#which input files missed chr7,8,X,Y?
Number of genes: 25735
Started run at 2025-10-02 22:25:46.767140
Finished at    2025-10-03 00:26, ~ 2 hours

##### checking total coverage for each DG:
awk '{sum=0; for (i=13; i<=42; i++) {sum+=$i;} print sum;}' data.txt


"""






################################################################################
def divbam(geneinfo,bam,sampleset):
    #process one gene, count reads for each sample from the input bam
    f = pysam.AlignmentFile(bam,'r'); DGdictdiv = {} #{DGID:{sample:count},...}
    #print(f)
    chrom,start,end,score,strand = geneinfo
    reads = f.fetch(chrom,start,end)
    for read in reads:
        QNAME,DGID = read.query_name,read.get_tag("DG") #print(QNAME,DGID)
        s = '-'.join(QNAME.split('-')[1:])
        if s not in sampleset: continue
        if DGID not in DGdictdiv: DGdictdiv[DGID] = {}
        if s not in DGdictdiv[DGID]: DGdictdiv[DGID][s] = 0
        DGdictdiv[DGID][s] += 1
    f.close()
    return DGdictdiv 
################################################################################





################################################################################
def exportDGs(geneinfo,DGdictdiv,DGbb,sampleset,outf):
    #DGdictdiv: {DGID:{sample:count},...}
    #directly export to another bed file for later conversion to bigwig. 
    chrom,start,end,score,strand = geneinfo
    f1 = pyBigWig.open(DGbb,'r'); genome = f1.chroms(); print(genome)
    if chrom not in genome: return 
    wigs = f1.entries(chrom,start,min(genome[chrom],end)) #DG bigbed, 1 per DG
    if not wigs: return
    for wig in wigs:
        s1,e1,string = wig; fields = string.split();
        DGID = ','.join(fields[0].split(',')[:3])
        if DGID not in DGdictdiv: continue
        str12 = '\t'.join([chrom,str(s1),str(e1),string])
        data = [DGdictdiv[DGID].get(s,0) for s in sampleset]
        strdata = '\t'.join([str(i) for i in data])
        outf.write(str12+'\t'+strdata+'\n')
    f1.close()
    return
################################################################################





################################################################################
if __name__ == "__main__":
    import os; os.environ['OPENBLAS_NUM_THREADS'] = '1' #override a bug on CARC
    import sys,pysam,pyBigWig
    from datetime import datetime
    from struclib import readbed6,makesamples

    if len(sys.argv) < 4: sys.exit("python bamDGdiv.py genesbed bam bed12")
    genes,bam,DGbb = sys.argv[1:4]; outf = open(bam[:-4]+'_DGdiv.bed','w')
    samples = makesamples(); sampleset = set(samples) #print(len(samples))
    genesdict = readbed6(genes)
    #for s in samples: print(s) #example: HepG2HNRNPCassomul

    print(str(datetime.today())+'\tProcessing all the genes')
    for gene in genesdict:
        print(str(datetime.today())+'\tProcessing', gene)
        geneinfo = genesdict[gene]
        DGdictdiv = divbam(geneinfo,bam,sampleset) #count reads for each DG
        exportDGs(geneinfo,DGdictdiv,DGbb,sampleset,outf) #export DG info to bed
    outf.close()
################################################################################



