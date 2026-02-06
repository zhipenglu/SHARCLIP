"""
kmer.py. zhipengluchina@gmail.com. 2025-06-17.
see clipmotifs.py on the selection of regions as peaks and backgrounds. 
input: two fasta files
output: count in a list and a volcano plot of ratio and counts. 

Worked very well. Showed clear motifs for HNRNPC in eCLIP data
The enrichment is much weaker for SHARCLIP, demonstrating reduced bias. 
Google: label points python scatter plot disperse
solutions to nicely label the dots on the plot. Also label all those with 4Ts. 

cd /Users/lu/Desktop/bedgraphs; pre=HepG2_HNRNPC_asso_indep_subtract_sorted
python ~/Dropbox/_scripts/kmer.py \
${pre}_all_background_1M.fa ${pre}_all_target_1M.fa 5 \
${pre}_volcano.txt ${pre}_volcano.pdf
#here background and target were from ...

cd /Users/lu/Desktop/eCLIP_desktop; pre=HepG2_HNRNPC_2_geonorm_all
python ~/Dropbox/_scripts/kmer.py \
${pre}_background_404374.fa ${pre}_target.fa 5 \
${pre}_5mer_volcano.txt ${pre}_5mer_volcano.pdf

pre=HepG2_RBFOX2_2_geonorm_all
python ~/Dropbox/_scripts/kmer.py \
${pre}_background_106848.fa ${pre}_target.fa 5 \
${pre}_5mer_volcano.txt ${pre}_5mer_volcano.pdf

pre=HepG2_TIA1_2_geonorm_all
python ~/Dropbox/_scripts/kmer.py \
${pre}_background_129436.fa ${pre}_target.fa 5 \
${pre}_5mer_volcano.txt ${pre}_5mer_volcano.pdf
"""




################################################################################
##### 1. function to count kmers in a sequence file.
def countk(file,k):
    #only count 4 capital letters "ACTG". Others are capitalized or ignored. 
    length,seq = 0,[]; list1 = ['ACTG' for i in range(k)] #
    kmers = {''.join(i):0 for i in list(product(*list1))} #kmer:count
    f = open(file, 'r')
    for line in f:
        if line[0] == '>':
            if seq:
                string = ''.join(seq).upper().replace('U','T'); seq = []
                for i in range(len(string)-k+1):
                    if string[i:i+k] in kmers: kmers[string[i:i+k]] += 1
        else: s = line.strip(); seq.append(s); length += len(s)
    f.close()
    string = ''.join(seq).upper().replace('U','T')
    for i in range(len(string)-k+1):
        if string[i:i+k] in kmers: kmers[string[i:i+k]] += 1
    return kmers, length #kmer dict {kmer:count}; length is sum of all seqs
################################################################################





################################################################################
##### 2. plot the kmers in a volcano plot.
def plotkmer(kmers1,len1,kmers2,len2):
    #A. plot all datapoints. 
    klist = list(kmers1.keys())
    freqs = [((kmers1[i]+1)/len1+(kmers2[i]+1)/len2)/2 for i in klist]
    logratios = [math.log2((kmers2[i]+1)/len2/(kmers1[i]+1)*len1)for i in klist]
    fig, ax = plt.subplots(); texts = []
    ax.scatter(logratios,freqs,color="b")
    
    #B. custom coloring and labeling for trimers. Add if needed. 
    labeldict = {"TTT":'y'; "TTTT":'r'}; texts = [] #yellow or red.
    for motif in labeldict: 
        ids1 = [i for i in range(len(klist)) if motif in klist[i]]
        klist1 = [klist[i] for i in ids1]
        freqs1,logratios1 = [freqs[i]for i in ids1], [logratios[i]for i in ids1]
        ax.scatter(logratios1,freqs1,color=labeldict[motif])
        for i,txt in enumerate(klist1):
            texts.append(ax.text(logratios1[i],freqs1[i],txt.replace("T","U")))
    adjust_text(texts,ax=ax); ax.set_xlim(-3.5,3.5); ax.set_ylim(0,0.005)
    plt.savefig(outpdf); plt.close()

    #C. export counts and frequencies to a text file
    outf = open(outtxt, 'w')
    for i in range(len(klist)):
        outf.write('\t'.join([klist[i],str(logratios[i]),str(freqs[i])])+'\n')
    outf.close()
    #for i in kmers: print(i,'\t',str(kmers[i]))
    #print(sum(list(kmers.values())),length,kmers)
    return 
################################################################################





################################################################################
if __name__ == "__main__":
    import sys, math
    from matplotlib import pyplot as plt
    from itertools import product
    from adjustText import adjust_text
    if len(sys.argv) < 6: sys.exit("python kmer.py fa1 fa2 k out.txt out.pdf")
    fa1,fa2 = sys.argv[1:3]; k = int(sys.argv[3]); outtxt,outpdf = sys.argv[4:6]
    kmers1,len1 = countk(fa1,k); kmers2,len2 = countk(fa2,k)
    plotkmer(kmers1,len1,kmers2,len2)
################################################################################





