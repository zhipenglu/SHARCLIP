"""
plotphaseRBP.py, zhipengluchina@gmail.com, 2025-07-17. 

Summary:
This analysis showed clear phasing of RBP binding sites on the HNRNPC arrays.
Next test the parameters for clustering, visualization, and build models. 
Some eCLIP files did not pass QC and were removed, e.g. HepG2/K562_RPS5_1/2.bam.

Each group of files: take ratio of 1 vs. 0 and 2 vs. 0. For example: 
HepG2_HNRNPC_asso_gap1_sorted_shortcntrs_K562_ZRANB2_0/1/2_xl.txt
Export individual normalized profiles and a composite heatmap ranked by peak pos

##### Caveats of the analyses listed above. 
1. Smoothening may reduce resolution, however, 5nt windows are small enough. 
2. Median normalization ignores the possibility that some proteins may bind
   a region very strongly or weakly.
    
##### 1. Example command on local computer:
mamba activate crssant
folder=/Users/lu/Desktop/eCLIP_desktop/RBPphase_K562masked/
pre=HepG2_HNRNPC_asso_gap1_sorted_shortpositioning_cntrs
pre=K562_HNRNPC_asso_gap1_sorted_shortpositioning_cntrs
python ~/Dropbox/_scripts/plotphaseRBP.py $folder $pre

"""


import sys, os
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib
from scipy import signal
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
#matplotlib.rcParams['pdf.fonttype'] = 42 
#plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['pdf.fonttype'] = 'none'


################################################################################
##### 1. group files 
def groupfiles(folder,cell):
    #example: HepG2_HNRNPC_asso_gap1_sorted_shortcntrs_K562_ZRANB2_0_xl.txt                                       100% 4153   423.1KB/s   00:00
    files = os.listdir(folder)
    fdict = {} #fdict: {(cell,RBP):{0:file0,1:file1,2:file2},...}
    for file in files:
        if prefix not in file: continue
        info = file.split('_'); cell,RBP,exp = info[-4:-1]
        if (cell,RBP) not in fdict: fdict[(cell,RBP)] = {}
        fdict[(cell,RBP)][int(exp)] = file
    #print(len(fdict.keys()))
    return fdict
################################################################################




################################################################################
##### 2. extract data from each file and normalize all samples against controls
def normalize(fdict,covmin,winsize,normmedian,smoothen):
    #fdict: {(cell,RBP):{0:file0,1:file1,2:file2},...}
    normall = {} # {(cell,RBP,exp):[],...} #exp=1/2, value: normalized list
    for k in fdict:
        group = fdict[k]
        if len(group) == 3:
            f = open(group[0],'r')
            ctrl = [float(i) for i in f.readlines()]; f.close()
            f = open(group[1],'r')
            exp1 = [float(i) for i in f.readlines()]; f.close()
            f = open(group[2],'r')
            exp2 = [float(i) for i in f.readlines()]; f.close()            
        r = range(1,len(ctrl))
        if ctrl[0]<covmin or exp1[0]<covmin: continue
        normall[k+(1,)] = [exp1[i]*ctrl[0]/ctrl[i]/exp1[0] for i in r]
        if ctrl[0]<covmin or exp2[0]<covmin: continue
        normall[k+(2,)] = [exp2[i]*ctrl[0]/ctrl[i]/exp2[0] for i in r]    
    #A. normalize each sample against its median.
    if normmedian: 
        for k in normall:  
            M = np.median(normall[k]); normall[k] = [i/M for i in normall[k]]
    #B. smoothen data, e.g., in 5nt windows:
    if smoothen: 
        window = np.ones(winsize)/winsize
        for k in normall: normall[k]=np.convolve(normall[k],window,mode='same')
    return normall # {(cell,RBP,exp):[],...} #exp=1/2, value: normalized list
################################################################################





#discarded. 
################################################################################
##### 3. check quality based on how rough the data tracks are. ...
#alternatively, remove low quality data where the RT counts too low. 
def lowquality(normall):
    #B. remove noisy data tracks: how do we define them?
    #for each pair of neighbor nts, compute (max-min)/average
    #normall: {(cell,RBP,exp):[],...} #exp=1/2, value: normalized list
    filtered = {}
    for k in normall:
        data = normall[k]
        r = range(len(data)-1)
        rough = [abs(data[i+1]-data[i])*2/(data[i+1]+data[i]) for i in r]
        rough = sum(rough)/len(rough)
        print('_'.join([str(i) for i in k]), '\t', str(rough))
    return 
################################################################################

        


    
################################################################################
##### 3. plot selected samples
def plotsample(normall,sample):
    #{(cell,RBP,exp):[],...} #exp=1/2, value: normalized list
    k = sample.split('-'); k = tuple(k[:2]+[int(k[2])])
    plt.plot(normall[k]); plt.xlim(0,200); plt.ylim(0.8,1.2)
    plt.savefig(sample+".pdf"); plt.close()
    return 
################################################################################





################################################################################
##### 3. plot selected samples, one RBP at a time
def plotRBP(normall,winsize,cell,RBP): #normall: {(cell,RBP,exp):[],...}
    reps = [normall[k] for k in normall if k[0]==cell and k[1]==RBP]
    colors = [(0.3,0.3,1-i*0.1,1) for i in range(len(reps))]
    xs = list(range(len(reps[0])))[winsize:-winsize]
    for i in range(len(reps)):
        plt.plot(xs,reps[i][winsize:-winsize],color=colors[i])
    plt.xlim(0,200); plt.ylim(0.8,1.2);
    plt.savefig(cell+'_'+RBP+".pdf"); plt.close()
    return 
################################################################################





#discarded. automatic clustering does not work well. 
################################################################################
##### 4. cluster and plot all samples in a heatmap.
def plotall(normall,cell):
    #normall: {(cell,RBP,exp):[],...} #exp=1/2, value: normalized list
    ratios = pd.DataFrame.from_dict(normall)
    fig, ax = plt.subplots()
    map1 = sns.clustermap(ratios,z_score=None,standard_scale=None,
                          metric="correlation",
                          row_cluster=False, col_cluster=True,
                          xticklabels=True, cmap="Greys", vmin=0.9,vmax=1.1)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=0.2)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=0.2)
    map1.savefig("RBPphase_heatmap.pdf"); plt.close()
    return 
################################################################################





################################################################################
##### 5. sort by nearest leftside peak center and plot the data:
def plotheatmap(normall,cell,winsize):
    #A. compute ranking statistic
    rankstats = []
    for k in normall: #normall: {(cell,RBP,exp):[],...}
        if k[0] != cell: continue
        d = normall[k]
        maxvd = float(max(d[winsize:-winsize])-min(d[winsize:-winsize])) #no use
        params = {"distance":30,"prominence":maxvd/3,"width":(5,40),"wlen":100}
        peaks,properties = signal.find_peaks(d,**params)
        ipsl,ipsr = properties["left_ips"],properties["right_ips"]
        proms = properties["prominences"]
        cs = [int((ipsl[i]+ipsr[i])/2) for i in range(len(ipsl))]
        left =  [c for c in cs if c <= 100]
        right = [c for c in cs if c > 100]
        index = min(right) if not left else max(left)
        index = max(left) if not right else min(right)
        rankstats.append(k+(maxvd,index))

    """
    #B. sort samples by the average indices of two samples for each RBP.
    #this ranking increased noise due to difference between replicates. Not used
    rankmerge = {} #key=RBP,value=[(cell,RBP,exp,maxvd,index),...]
    RBPs = set([i[1] for i in rankstats])
    for info in rankstats: #[((cell,RBP,exp),maxvd,index),...]
        if info[1] not in rankmerge: rankmerge[info[1]] = []
        rankmerge[info[1]].append(info) #samples of one RBP merged together
    ranklist = list(rankmerge.values()) #print(ranklist)
    #below: x is (cell,RBP,exp,maxvd,index)
    rank = sorted(ranklist,key=lambda x:sum([i[4] for i in x])/len(x))
    rank = [j for i in rank for j in i] #flatten list. 
    """
    #B. sort based on peak positions closest to the center of the plotting area.
    rank = sorted(rankstats,key=lambda x:x[4]) 
    
    #C. plot the heatmap
    normall1 = {k:normall[k[0:3]] for k in rank}
    data = pd.DataFrame.from_dict(normall1)
    fig, ax = plt.subplots()
    map1 = sns.heatmap(data,cmap="Greys", vmin=0.9,vmax=1.1);
    fig = map1.get_figure(); fig.savefig(cell+"_RBPphase_heatmap.pdf");
    plt.close()
    print("Plotted samples:", len(normall))
    return
################################################################################





################################################################################
##### 5. process all files
if __name__ == "__main__":
    #A. group eCLIP files
    if len(sys.argv)<3: sys.exit("python plotphaseRBP.py folder prefix")
    folder,prefix = sys.argv[1:3]; cell = prefix.split('_')[0]
    fdict = groupfiles(folder,prefix) #{(cell,RBP):{0:file0,1:file1,2:file2}}
    fdict = {k:fdict[k] for k in fdict if len(fdict[k])==3} #remove bad samples
        
    #B. extract, normalize and smoothen the data for each group. 
    covmin = 1E6 #min RT count in a sample to disquality the dataset
    winsize = 5  #for smoothening the ratio curves. 
    normall = normalize(fdict,covmin,winsize,normmedian=True,smoothen=True)

    #C. plot heatmap after sorting by first peak to the left of figure center.
    plotheatmap(normall,cell,winsize)  

    #D. plot each RBP for a specific cell line
    #plotsample(normall,sample) #individual samples
    RBPs = set([k[1] for k in normall if k[0]==cell])
    for RBP in RBPs: plotRBP(normall,winsize,cell,RBP)
    
################################################################################








