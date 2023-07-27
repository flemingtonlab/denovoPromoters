#!/usr/bin/env python

import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
import random

# Look for ATAC coverage over Zta CHIP peaks (RAJI cells)
# Plots created 5/4/2023


# Zta ChIP 
bed_path = '/Users/nate/Documents/Projects/De_Novo_Promoter/de_novo_promoter_heatmap/Raji_Zta_ChIP_induced_pooledReps_summits.hg38.bed'

xx = pd.read_table('/Users/nate/Raji_Zta_ChIP_induced_summits_slop5000.no_blacklist.bed', index_col=0, header=None)
xx = set(xx[3])


akata_uninduced_str1_paths = [
    '/Users/nate/dnovo/atac/akata/deduplicated_Akata-C1_hg38chr1-22XYplus_chrEBVAkataInverted.bw',
    '/Users/nate/dnovo/atac/akata/deduplicated_Akata-C2_hg38chr1-22XYplus_chrEBVAkataInverted.bw',
    '/Users/nate/dnovo/atac/akata/deduplicated_Akata-C3_hg38chr1-22XYplus_chrEBVAkataInverted.bw',
    ]

akata_induced_str1_paths = [
    '/Users/nate/dnovo/atac/akata/deduplicated_Alata-BCR1_hg38chr1-22XYplus_chrEBVAkataInverted.bw',
    '/Users/nate/dnovo/atac/akata/deduplicated_Alata-BCR2_hg38chr1-22XYplus_chrEBVAkataInverted.bw',
    '/Users/nate/dnovo/atac/akata/deduplicated_Alata-BCR3_hg38chr1-22XYplus_chrEBVAkataInverted.bw'
    ]

mutu_ctl_str1_paths = [
    '/Users/nate/dnovo/atac/mutu/deduped_MC2_hg38AkataInverted.filtered.bw',
    '/Users/nate/dnovo/atac/mutu/deduped_MC4_hg38AkataInverted.filtered.bw',
    '/Users/nate/dnovo/atac/mutu/deduped_MC5_hg38AkataInverted.filtered.bw'
]


mutu_zta_str1_paths = [
'/Users/nate/dnovo/atac/mutu/deduped_MZ1_hg38AkataInverted.filtered.bw',
'/Users/nate/dnovo/atac/mutu/deduped_MZ2_hg38AkataInverted.filtered.bw',
'/Users/nate/dnovo/atac/mutu/deduped_MZ4_hg38AkataInverted.filtered.bw'
]

akata_uninduced_total_coverages = [76183920, 74468502, 58841542]
akata_induced_total_coverages = [22406826, 18556596, 19483240]
mutu_ctl_total_coverages = [67952942, 53139584, 53456816]
mutu_zta_total_coverages = [70348786, 65719606, 69502884]

stranded = False

upstream_bases = 5000
downstream_bases = 5000
length = upstream_bases + downstream_bases

regions = []
with open(bed_path) as bed_handle:
    for line in bed_handle:
        regions.append(line.strip('\n').split('\t'))

random.shuffle(regions)
regions.sort(key=lambda x:float(x[4]))
regions = [i for i in regions if i[3] in xx]


def extract_coverage(strand1_paths, total_coverages, regions):
    
    for str1_path, total_coverage in zip(strand1_paths, total_coverages):
        str1 = pyBigWig.open(str1_path)
        m = np.zeros([len(regions), length])

        for ind, region in enumerate(regions):
            start = int(float(region[1]))
            stop = int(float(region[2])) 
            #strand = region[5]

            middle = (start + stop) // 2 
            start = middle - upstream_bases
            stop = middle + downstream_bases
            vals = str1.values(region[0], start, stop)
            
            vals = np.abs(vals)
            vals = np.nan_to_num(np.array(vals), nan=.01)

            #if strand == '-':
             #   vals = np.flip(vals)
            vals += .01

            m[ind] += (10000000 * vals)/ total_coverage
        
        print(str1_path, "done")
    return m / len(strand1_paths)

def get_row_cov_max(matrices_list, min_divisor=1):
    rows = matrices_list[0].shape[0]
    number_of_matrices = len(matrices_list)
    maxes = np.zeros([rows, number_of_matrices + 1])
    for col, matr in enumerate(matrices_list):
        maxes[:, col] = np.max(matr, 1)
    maxes[number_of_matrices] = min_divisor
    return np.max(maxes, 1)


    
akata_uninduced_coverage = extract_coverage(akata_uninduced_str1_paths, akata_uninduced_total_coverages, regions=regions)
akata_induced_coverage = extract_coverage(akata_induced_str1_paths, akata_induced_total_coverages, regions=regions)
mutu_ctl_coverage = extract_coverage(mutu_ctl_str1_paths, mutu_ctl_total_coverages, regions=regions)
mutu_zta_coverage = extract_coverage(mutu_zta_str1_paths, mutu_zta_total_coverages, regions=regions)

height = 8 
all_matrices = [akata_uninduced_coverage, akata_induced_coverage, mutu_ctl_coverage, mutu_zta_coverage]
for i,j in zip(all_matrices, ["ak_un", "ak_in", "mut_ctl", "mut_zta"]):
    fig = plt.figure(figsize=(3, height))
    plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian", vmax=1)
    plt.xticks([length//2])
    plt.yticks([])
    #plt.axvline(length//2, ls='--',c='k')
    plt.savefig(j + f'_denovo_noblacklist_atac_{length}.zta_raji_chip_sites.svg')
    plt.close('all')

fig = plt.figure()
akun_normed = np.sum(akata_uninduced_coverage, 0) / len(regions)
akin_normed = np.sum(akata_induced_coverage, 0) / len(regions)
plt.plot(range(length), akun_normed)
plt.fill_between(range(length), akun_normed, alpha=.3)
plt.plot(range(length), akin_normed)
plt.fill_between(range(length), akin_normed, alpha=.3)    
plt.xticks([length//2])
plt.yticks([])
plt.xlim([0,length])
#plt.ylim([0, 1.2])
#plt.axvline(length//2, ls='--',c='k')
plt.savefig(f'akata_summation_curve_atac_noblacklist_{length}.zta_raji_chip_sites.svg')
plt.close('all')

fig = plt.figure()
muctl_normed = np.sum(mutu_ctl_coverage, 0) / len(regions)
muzta_normed = np.sum(mutu_zta_coverage, 0) / len(regions)
plt.plot(range(length), muctl_normed)
plt.fill_between(range(length), muctl_normed, alpha=.3)
plt.plot(range(length), muzta_normed)
plt.fill_between(range(length), muzta_normed, alpha=.3)    
plt.xticks([length//2])
plt.xlim([0,length])
#plt.ylim([0, 0.8])
plt.yticks([])
#plt.axvline(length//2, ls='--',c='k')
plt.savefig(f'mutu_summation_curve_atac_noblacklist_{length}.zta_raji_chip_sites.svg')
plt.close('all')