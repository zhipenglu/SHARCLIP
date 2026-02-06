"""
motifperiodicity.py, zhipengluchina@gmail.com, 2025-07-06. 
Check intermotif distance distribution for (N)m of A, C, G and T (not U), m=3-10
The first nucleotide position on the RNA strand is used to calculate distance.
Only consider noncontinuous motifs, i.e. only the first occurrence is
counted in a long homopolymer. Analyzes one kmer at a time.
See script kmerdist.py for systematic analysis of all kmers k=4/5. 
Only processed the genes on the 25 normal chromosomes, 1-22, X, Y and M
Shuffling test is better than theoretical distribution but time consuming. 

##### python histogram implementations: 
#Left-inclusive, right-exclusive: [a, b), a <= x < b 
#Last bin exception: [c, d], c <= x <= d
##### sum(list,[]) is very slow, O(L**2). See: https://stackoverflow.com/\
questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
Instead, use this or other methods: [item for sublist in l for item in sublist]

##### Conclusions:
1. There is a clear spike at ~154nt, consistent with the ~39nt periodicity
1. The 39nt periodicity may remain after removing the Alu elements. 
2. This 154nt periodicity is shorter than Konig 2011. 
3. Separation of intron vs. exons, ... Most exons are not long enough
4. The longer the T stretch, the more obvious the ~150-160nt and other peaks.
5. Is it due to the asymmetry of the tetramer? May be C2 is needed for binding?
6. The ~290nt peak is likely linked to the ~135nt peak, with a distance of ~154
If the 154/287nt periodicity is caused by the millions of Alus, we are largely
overturning the Konig model. The Alus cannot force HNRNPC tetramers to evolve a
39nt periodicity. Exclude the Alu tetramers and then check for motifs.


################################################################################
##### all 24 chromosomes separately. ONLY FOR DEBUGGING. 
for i in $(seq 1 22); \
do (grep -w chr$i hg38refGene.bed > hg38refGene_chr$i.bed &); done
grep chrX hg38refGene.bed > hg38refGene_chrX.bed
grep chrY hg38refGene.bed > hg38refGene_chrY.bed
script=~/Dropbox/_scripts/motifperiodicity.py
f1=~/Documents/lulab/projects/hg38/fasta
f2=~/Documents/lulab/projects/hg38/refGene/hg38refGene
for i in $(seq 1 22); \
do (python $script ${f1}/chr$i.fa ${f2}_chr$i.bed TTTT chr$i &); done
python $script ${f1}/chrX.fa ${f2}_chrX.bed TTTT chrX
python $script ${f1}/chrY.fa ${f2}_chrY.bed TTTT chrY


################################################################################
##### entire chr11
genes=~/Documents/lulab/projects/hg38/refGene/hg38refGene_chr11.bed
#mask=no
fasta=~/Documents/lulab/projects/hg38/fasta/chr11.fa 
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes yes no 4 4 800 &
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes no no 4 4 800 &
#mask=yes, repeats in lower case
fasta=/Users/lu/igv/genomes/hg38/hg38.fa 
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes yes yes 4 4 800 &
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes no yes 4 4 800 &


################################################################################
##### entire genome, only 25 normal chromosomes. 
genes=~/Documents/lulab/projects/hg38/refGene/hg38refGene_25chroms.bed
#mask=no
fasta=~/Documents/lulab/projects/hg38/fasta/hg38_25chroms.fa 
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes yes no 3 10 800 &
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes no no 3 10 800 &
#mask=yes, repeats in lower case 
fasta=/Users/lu/igv/genomes/hg38/hg38.fa 
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes yes yes 3 10 800 &
python ~/Dropbox/_scripts/motifperiodicity.py $fasta $genes no yes 3 10 800 &


"""



import sys, os, itertools, random, datetime
from itertools import product
from matplotlib import pyplot as plt
from struclib import rc,readfa,readbed6
if len(sys.argv)<5:
    sys.exit("python motifperiodicity.py fasta genes shuffle mask \
             lmin lmax plotlen")
    #fasta, repeats in lower or upper case.
    #genes in bed6 format, introns or exons. 
    #shuffle: yes or no.
    #soft mask: yes or no, masked means repeats in fasta is lower case
    #lmin,lmax, motif length min and max, typically 3-10. 
    #plotlen, default 800. 
#readfa(fa), returns fadict
#readbed6(bed6), returns genesdict: {name:[chrom,start,end,score,strand], ...}.




################################################################################
##### 1. make the following motifs
def makemotifs(lmin,lmax):
    homo = {nt*i:0 for i in range(lmin,lmax+1) for nt in "ATCG"}   #32
    #kmer4 = {''.join(i):[] for i in product("ATCG", repeat=4)}  #256
    #kmer5 = {''.join(i):[] for i in product("ATCG", repeat=5)}  #1024
    return homo
################################################################################





################################################################################
##### 2. find motif positions in the genome
def findmotifs(fadict,genesdict,motif,chroms,maxlen,flag):
    """
    fadict/genesdict: {chrom:seq,...},{name:[chrom,start,end,score,strand],...}
    motif: e.g. TTTT; chroms: list of chromosomes; maxlen=1E9 for now
    flag: shuffle or not (yes vs. no)
    """
    allpos = []; counter = 0 #allpos: list of lists of positions for each chrom
    for gene in genesdict:
        counter += 1
        if not counter%1000:
            time = str(datetime.datetime.today())
            print(time + " Processed {} genes".format(counter))
        chrom,start,end,score,strand = genesdict[gene]
        loci = []; mlen = len(motif); glen = end-start
        if chrom not in fadict or chrom not in chroms or glen>maxlen: continue
        if strand == '-': motif = rc(motif)
        seq = fadict[chrom][start:end]
        if flag: random.seed(0); l=list(seq); random.shuffle(l); seq=''.join(l)
        for i in range(len(seq)-mlen):
            if seq[i:i+mlen]==motif: loci.append(i if strand=="+" else i+mlen-1)
        allpos.append(loci)
    print("Total number of processed genes:", counter)
    return allpos
################################################################################





################################################################################
##### 3. export inter-motif distances to txt and plot in pdf
def exportdist(genesbed,allpos,motif,plotlen,flag,mask):
    dists = []
    for x in allpos: dists.append([x[i+1]-x[i] for i in range(len(x)-1)])
    dists = [i for dist in dists for i in dist if i>1]
    fig, ax = plt.subplots(); ranges = list(range(plotlen))
    n, bins, patches = plt.hist(dists,ranges,density=True,facecolor="blue")
    figname = genesbed.split('/')[-1][:-4] + "_motif_{}.pdf".format(motif)
    if flag: figname = figname[:-4] + "_shuffle.pdf"
    if mask: figname = figname[:-4] + "_mask.pdf"
    plt.savefig(figname); plt.close()
    outf = open(figname[:-4]+"_dist.txt",'w')
    for i in range(plotlen-1): outf.write(str(n[i])+'\n')
    print(str(datetime.datetime.today())+"\tFinished motif {}".format(motif))
    outf.close()
    print("Number of inter-motif distances:",len(dists))
    return
################################################################################






################################################################################
##### 4. plot individual motifs for each collection of motifs
if __name__ == '__main__':
    fasta,genesbed,shuffle,mask,lmin,lmax,plotlen = sys.argv[1:8]
    chroms = ["chr"+str(i) for i in range(1,23)] + ["chrX","chrY"]
    maxlen = 1E9 #ignore superlong genes. Use 1E9 to remove limit. 
    plotlen = int(plotlen) #longer than a monoparticle
    flag = 1 if shuffle == "yes" else 0
    mask = 1 if mask == "yes" else 0
    fadict = readfa(fasta); genesdict = readbed6(genesbed)
    homo = makemotifs(int(lmin),int(lmax)) #only analyzes the homo polymers
    for motif in homo: 
        allpos = findmotifs(fadict,genesdict,motif,chroms,maxlen,flag)
        dists = exportdist(genesbed,allpos,motif,plotlen,flag,mask)
################################################################################




