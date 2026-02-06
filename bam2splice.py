"""
sam2splice.py, zhipengluchina@gmail.com, Written 2022-05-03. Updated 2024-01-29.
rewritten on 2025-08-08 to speed up processing, using unique counts. 
Note: f346N indicates fields of 3, 4 and 6 for records with N-containing CIGARs. 
A simple script to visualize splice junctions. The arcs in the bed file can be
colored/gray-scaled based on junction read count. This approach is cleaner than
Sashimi. See R01 grant Figs. 7/8/12. Manually label junction counts for now.
Alternatively, make another bed without arcs, where name is the count. 
In rMATS, --anchorLength is set to 1 as the default. We use the same here. 

Useful files:
1. original bam files
2. f346 converted bed files for arcs
3. coverage bw files.
4. rMATS output statistics. Made volcano figures already. p<0.01, dpsi>0.2

Input1: 3 fields from BAM: RNAME POS CIGAR, where CIGAR contains N
Input2: junction database to filter splicing events. Now using RJunBase. 
Output: splice junctions in an arc bed format.

##### 1. Extract the 3 fields from BAM: original analysis on ENCODE data.
#e.g. HepG2_hg38all_HNRNPH1, 269506 reduced to 15966 lines.
#did not work well on transcriptome due to singletons, e.g. 856M to 710M
bam=HepG2/K562_hg38all.bam;
f346=${bam/.bam/_f346.txt}; f346c=${f346/.txt/_count.txt}
salloc samtools view $bam | awk '$6~/N/ {print $3 "\t" $4 "\t" $6}' >$f346 &
salloc sort $f346 | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' >$f346c &
HepG2_hg38all.bam: 76G; K562_hg38all.bam: 82G
Our datasets: 5 cell lines, 4 replicates each, ~5G bam, no need for splitting.

##### 2. Use a dict e.g. GENCODE or RJuncBase to filter junctions, Li 2021 NAR.
682018 linear, back and fusion junctions from ~30k normal/cancerous samples. 
/Users/lu/Documents/lulab/projects/hg38/RJunBase/detail_LS_annotation.txt

##### 3. bam2splice.py, parse CIGAR strings, build splicedict and write output
cd ~/Documents/lulab/projects/CLIP/encode/encodesplice
script=~/Dropbox/_scripts/bam2splice.py
f346c=HepG2_hg38all_HNRNPH1_count.txt; bedarc=${f346/txt/bed}
juncdb=~/Documents/lulab/projects/hg38/RJunBase/detail_LS_annotation_HNRNPH1.txt
python $script $f346c $juncdb $bedarc
#this run results in 42 detected splice junctions in HepG2_hg38all.

##### 4. example command on CARC:
cd /project/zhipengl_72/zhipengl/encode
script=/project/zhipengl_72/labshare/scripts/bam2splice.py
bam=HepG2/K562_hg38all.bam
f346=${bam/.bam/_f346.txt}; bedarc=${f346/.txt/_splice.bed}
juncdb=/project/zhipengl_72/labshare/hg38/detail_LS_annotation.txt
salloc --mem=10G --time=10:00:00 python $script $f346 $juncdb $bedarc &

"""





################################################################################
#####1. Build a splice dictionary from GENCODE or RJunBase #how big is it 
def makesplicedict(juncdb): 
    juncf = open(juncdb, 'r')
    juncdict = {} # (RNAME,start,end): gene_name; strand ignored here. 
    for line in juncf:
        record = line.split()
        if record[0] == "JunctionID": continue
        gene, junc = record[0].split('_')[-2], record[1]
        RNAME, coords, strand = junc.split(':')
        coords = [int(i) for i in coords.split('|')]
        juncdict[(RNAME,coords[0],coords[1])] = gene
    juncf.close()
    return juncdict
################################################################################





################################################################################
#####2. Function: convert each input read info to a list of splices
def makesplices(RNAME, POS, CIGAR): 
    #gaps, a list of N strings, e.g. ['100N', '200N']
    #segs, a list of segments, each with all possible operations except N
    #gaplens, a list of N lengths, e.g. [100, 200]
    #Rlens, segment lengths of consumed Reference [MD=X]
    gaps=re.findall(r'\d+N', CIGAR)
    segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
    gaplens=[int(gap[:-1]) for gap in gaps] #gap lengths 
    Rlens=[sum([int(i[:-1]) for i in re.findall(r'\d+[MD=X]',s)]) for s in segs]
    splices = [] #[[RNAME, LEFT, RIGHT, LEN], ...]
    for i in range(len(gaps)): #n gaps and n+1 segs
        splices.append((RNAME, POS+sum(gaplens[:i])+sum(Rlens[:i+1])-1, \
                       POS+sum(gaplens[:i+1])+sum(Rlens[:i+1])))
    return splices
################################################################################





################################################################################
#####3. Extract and count splice junctions, export to bedarc file
def countsplices(f346,juncdict): 
    f346f = open(f346,'r'); i = 0
    splicedict = {} #{(RNAME,pos1,pos2):[gene,count,...}
    for line in f346f:
        i += 1
        if not i%1000000: print("Processed", i, "lines ...")
        RNAME,POS,CIGAR = line.split(); POS = int(POS)
        if "N" not in CIGAR: continue #in case non-filtered
        splices = makesplices(RNAME,POS,CIGAR)
        for splice in splices:
            if splice in juncdict:
                if splice not in splicedict:
                    splicedict[splice] = [juncdict[splice],1]
                else: splicedict[splice][1] += 1
    f346f.close()
    return splicedict
################################################################################





################################################################################
#####1. write splicedict to bedarc: a bed file to display arcs. 
def writebed(splicedict,bedarc,header): 
    bedarcf = open(bedarc, 'w'); bedarcf.write(header)
    #A. rearrange splicedict to genedict, where key=gene, for normalization
    genedict = {} #{gene: [[gene,(RNAME,start,end),count], ...],...}
    for splice in splicedict: #{(RNAME,start,end):[gene,count],...} 
        gene,count = splicedict[splice] #splice: (RNAME,start,end)
        if gene not in genedict: genedict[gene] = []
        genedict[gene].append([gene,splice,count])
        
    #B. add normalized RGB color and export to bed arc format
    for gene in genedict: 
        splices = genedict[gene]; Max = max([i[2] for i in splices])
        for i in range(len(splices)):
            gene,splice,count = splices[i]
            RGB = ','.join([str(int(230*(1-count/Max)))]*3)
            record = list(splice)+[count,gene,'.',splice[1],splice[2],RGB]
            bedarcf.write('\t'.join([str(i) for i in record]) +'\n')
    bedarcf.close()
    return 
################################################################################





################################################################################
if __name__ == "__main__": 
    import sys, re
    from datetime import datetime
    if len(sys.argv)<4:
        print("Usage: python bam2splice.py f346 juncdb bedarc")
        print("f346: 3 fields from BAM: RNAME POS CIGAR")
        print("bedarc: bed file to visualize the arcs")
        print("requires the starting BAM file to be sorted")
        sys.exit()
    f346,juncdb,bedarc = sys.argv[1:4]; header = "track graphType=arc\n"
    print(str(datetime.today())+'\tReading splice database')
    juncdict = makesplicedict(juncdb)        #1. get junction annotations
    print(str(datetime.today())+'\tCounting splice junctions')
    splicedict = countsplices(f346,juncdict) #2. count splicing events
    print(str(datetime.today())+'\tWriting output bed file')
    writebed(splicedict,bedarc,header)       #3. write to bedarc file format
    print(str(datetime.today())+'\tFinished analysis')
    #print("Junction DB size:", len(makesplicedict(juncdb))) #n=682018
################################################################################








