
import sys
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import random
import re
import glob

bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.noblacklist.bed'
canonical_bed_path = '/Users/nate/dnovo/beds/TSS_from_mutu_ctl2.bed'

akata_uninduced_str1_paths = ['/Volumes/de_novo/ATAC/Akata_ATAC/bams/AkataC24h.footprints.bw']
akata_induced_str1_paths = ['/Volumes/de_novo/ATAC/Akata_ATAC/bams/AkataIg24h.footprints.bw']

mutu_ctl_str1_paths = ['/Volumes/de_novo/ATAC/Mutu_ATAC/bams/MC24.footprints.bw']
mutu_zta_str1_paths = ['/Volumes/de_novo/ATAC/Mutu_ATAC/bams/MZ24.footprints.bw']

upstream_bases = 5000
downstream_bases = 5000
length = upstream_bases + downstream_bases


def extract_coverage(strand1_paths, regions):
    
    for str1_path in strand1_paths:
        str1 = pyBigWig.open(str1_path)
        m = np.zeros([len(regions), length])
        for ind, region in enumerate(regions):
            start = int(float(region[1]))
            stop = int(float(region[2])) 
            strand = region[5]
            middle = (start + stop) // 2 
            start = middle - upstream_bases
            stop = middle + downstream_bases
            vals = str1.values(region[0], start, stop)
            vals = np.abs(vals)
            vals = np.nan_to_num(np.array(vals), nan=0)
            if strand == '-':
                vals = np.flip(vals)
            m[ind] +=  vals
        print(str1_path, "done")
    return m / len(strand1_paths)

regions = []
with open(canonical_bed_path) as bed_handle:
    for line in bed_handle:
        regions.append(line.strip('\n').split('\t'))

random.shuffle(regions)
regions.sort(key=lambda x:float(x[6]))

        
mutu_ctl_coverage = extract_coverage(mutu_ctl_str1_paths, regions=regions)
mutu_zta_coverage = extract_coverage(mutu_zta_str1_paths, regions=regions)

height = 8 
all_matrices = [mutu_ctl_coverage, mutu_zta_coverage]
for i,j in zip(all_matrices, ["mut_c", "mut_z"]):
    fig = plt.figure(figsize=(3, height))
    plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian", vmax=.5)
    plt.xticks([length//2])
    plt.yticks([])
    plt.savefig(j + f'_denovo_noblacklist_atac_mutu_{length}.canonical.svg')
    plt.close('all')

# bedtools closest -a de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed -b CANONICAL_TSS_from_mutu_ctl.bed -s -id -D a>upoG_all_canonical.bed
# bedtools closest -a de_novo_promoters_Akata_BCR_plus_Mutu_Zta.sorted.bed -b CANONICAL_TSS_from_mutu_ctl.bed -S -iu -D a>antisense_allcanonical.bed
upog = pd.read_table('~/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure6/starting_files/upoG_all_canonical.bed', header=None)
anti = pd.read_table('~/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure6/starting_files/antisense_allcanonical.bed', header=None)

# Within 1.5kb of canonical start
upog = upog[(upog[12]>-1500) & (upog[12]<-1)]
anti = anti[(anti[12]<1500) & (anti[12]>-1)]

# High CAGE counts only (30)
upog30 = upog[upog[4]>30]
anti30 = anti[anti[4]>30]

# Genes nearby
both = set(upog30[10]) | set(anti30[10])

l = []

for i in regions:
    if i[4] in both:
        l.append('k')
    else:
        l.append('w')

fig = plt.figure(figsize=(.5,8))
ax = plt.subplot()
for i in range(len(l)):
    if l[i] == 'k':
        plt.plot([0,1],[i,i],c=l[i],lw=.1)

plt.ylim([0,len(l)])
plt.xlim([0,1])
plt.xticks([])
plt.yticks([])
plt.savefig('barcode.svg',dpi=300)


fig = plt.figure(figsize=(8,1))
ax = plt.subplot()
expression = np.array([float(i[6]) for i in regions])
y_vals = range(len(expression))
plt.xlim([0,len(expression)])
plt.ylim([np.min(np.log2(expression+.1)),np.max(np.log2(expression+.1))])
plt.xticks([])
plt.yticks([])
plt.fill_between(y_vals,np.log2(expression+.1),np.min(np.log2(expression+.1)), color='k')
plt.plot(y_vals, np.log2(expression+.1),c='k')
plt.savefig('expression.png')