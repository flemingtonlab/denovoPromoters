
import sys
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import random
import glob
import re

# Plot created May 17, 2023
# Mutu 12h ATAC samples
# Any de novo site that had overlap with a blacklisted region (encode) was removed (using the -5000 to +5000 region)

bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.noblacklist.bed'
H3K4ME3_str1_paths = ['/Volumes/de_novo/ChIP/Akata/bigwigs/H3K4M3Lat.merged.bam.bw']
H3K27ME3_str1_paths = ['/Volumes/de_novo/ChIP/Akata/bigwigs/H3K27M3Lat_1.head.fq.chrMT_chrEBV_unmapped_removed.bam.bw']
H3K27AC_str1_paths = ['/Volumes/de_novo/ChIP/Akata/bigwigs/H3K27AcLat_1.fq.gz.head.fq.chrMT_chrEBV_unmapped_removed.bam.bw']
H3K4ME3_REAC_str1_paths = ['/Volumes/de_novo/ChIP/Akata/bigwigs/H3K4M3Reac_1.merged.bam.bw']


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
            try:
                strand = region[5]
            except:
                strand = '+'
                
            middle = (start + stop) // 2 
            start = middle - upstream_bases
            stop = middle + downstream_bases
            vals = str1.values(region[0], start, stop)
            vals = np.abs(vals)
            vals = np.nan_to_num(np.array(vals), nan=0)
            if strand == '-':
                vals = np.flip(vals)
            
            m[ind] += vals  
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


def import_bed(path, col_sort=4):
    regions = []
    with open(path) as bed_handle:
        for line in bed_handle:
            regions.append(line.strip('\n').split('\t'))
    random.shuffle(regions)
    regions.sort(key=lambda x: float(x[col_sort]))
    return regions

def plot_heatmap(matrix, name, colormap='Reds',height=8):
    fig = plt.figure(figsize=(3,height))
    plt.imshow(matrix, cmap=colormap,aspect='auto', interpolation="gaussian", vmax=1)
    plt.xticks([])
    plt.yticks([])
    plt.savefig(name + '.heatmap.svg')
    plt.close('all') 

def plot_sumcurves(matrix1, name, color='r',cutoff_number=.0765, binsize=10, max_y=3.5):
    matrix1[matrix1 < cutoff_number] = 0
    matrix1[matrix1 > cutoff_number] = 1
    sums1 = np.sum(matrix1, 0)
    sum_hist1 = [0]
    for i in range(0, len(sums1), binsize):
        sum_hist1.append(np.sum(sums1[i:i+binsize]))
    sum_hist1.append(0)
    fig = plt.figure(figsize=(4,4)) 
    y_vals= np.array(sum_hist1[1:-1]) / matrix1.shape[0]
    plt.plot(y_vals, color=color)
    plt.fill_between(range(len(y_vals)), y_vals, y2=[np.min(y_vals)]*len(y_vals), color=color)
    plt.ylim([np.min(y_vals), max_y])
    plt.yticks([])
    plt.xticks([])
    print(np.max(sum_hist1[1:-1]) / np.shape(matrix1)[0],name)
    plt.savefig(name + '.summation.svg')

def get_zta_sites(denovo_regions):
    zta_vpic_df = pd.read_table('/Users/nate/Downloads/motifs_in_de_novo_promoters/4_vPIC_plus_Zta/4_TATTAAA_TATTTAA_and_Zta_ChIP_seq_Hammerschmidt/TATTAAA_plus_TATTTAA_and_Zta_ChIP_Hammerschmidt_output/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed_TATTAAA_plus_TATTTAA.bed.binary_motif_info.bed_Raji_Zta_ChIP_induced_pooledReps_summits.bed.no_EBV_plus_strand.bed.binary_motif_info.tsv', index_col=0)
    zta_vpic_df = zta_vpic_df.set_index('DN prom name')
    zta = zta_vpic_df[zta_vpic_df['Zta_ChIP_Hammerschmidt'] == 1]
    zta = set(zta.index)
    return [i for i in denovo_regions if i[3] in zta]

def get_vpic_sites(denovo_regions):

    bcrf1_motif = "TATT[TA]AA"
    bcrf1_prog = re.compile(bcrf1_motif)
    working_dir = "/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure3/starting_files/"
    fa40_25 = {}
    with open(working_dir + 'de_novo_tss_40to25bp_upstream.allsites.fa') as infile:
        for line in infile:
            line = line.strip('\n')
            if line[0] == '>':
                keya = line.replace('>','').split('(')[0]
            else:
                fa40_25[keya]= line
    return [i for i in denovo_regions if bcrf1_prog.search(fa40_25[i[3]]) is not None]
    
def get_rta_sites(denovo_regions):
        
    fa200_0 = {}
    working_dir = "/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure3/starting_files/"

    with open(working_dir + 'de_novo_tss_200bp_upstream.allsites.fa') as infile:
        for line in infile:
            line = line.strip('\n')
            if line[0] == '>':
                keya = line.replace('>','').split('(')[0]
            else:
                fa200_0[keya]= line
    rta_motif_fw = "G[ACTG]CC[ACGT]{8,10}GG[ACGT]G" 
    rta_motif_rev= "C[ACGT]CC[ACGT]{8,10}GG[ACGT]C"
    rta_prog = re.compile(f'{rta_motif_fw}|{rta_motif_rev}')
    return [i for i in denovo_regions if rta_prog.search(fa200_0[i[3]]) is not None]


def get_sites(denovo_regions):
    rta = get_rta_sites(denovo_regions)
    zta = get_zta_sites(denovo_regions)
    vpic = get_vpic_sites(denovo_regions)
    all_three = rta+zta+vpic
    all_three = [i[3] for i in all_three]

    l = []
    two_plus_sites = []
    for i in all_three:
        if i[3] in l:
            two_plus_sites.append(i[3])
        l.append(i[3])

    no_site = [i for i in denovo_regions if i[3] not in all_three]
    vpic = [i for i in vpic if i[3] not in two_plus_sites]
    zta = [i for i in zta if i[3] not in two_plus_sites]
    rta = [i for i in rta if i[3] not in two_plus_sites]
    print(len(rta) + len(zta) + len(vpic) + len(no_site) + len(set(two_plus_sites)), len(denovo_regions))

    return rta, zta, vpic, no_site

denovo_regions = import_bed(bed_path)
rta, zta, vpic, no_site = get_sites(denovo_regions)


for reg, reg_name in zip([rta, zta, vpic, no_site], ["rta", "zta", "vpic", "none"]):
            
    h3k4me3_coverage = extract_coverage(H3K4ME3_str1_paths, regions=reg)
    # h3k27me3_coverage = extract_coverage(H3K27ME3_str1_paths, regions=reg)
    # h3k27ac_coverage = extract_coverage(H3K27AC_str1_paths, regions=reg)
    h3k4me3_reac_coverage = extract_coverage(H3K4ME3_REAC_str1_paths, regions=reg)
    height = 8 * (len(reg) / len(denovo_regions))
    all_matrices = [h3k4me3_coverage, h3k4me3_reac_coverage]
    for i,j in zip(all_matrices, ["h3k4me3", 'h3k4me3_react']):
       
        fig = plt.figure(figsize=(3, height))
        plt.imshow(i, cmap='Purples',aspect='auto', interpolation="gaussian", vmax=100)
        plt.xticks([length//2])
        plt.yticks([])
        plt.savefig(j + f'_denovo_histonecov_{length}.{reg_name}.vmax100.svg')
        plt.close('all')


for reg, reg_name in zip([rta, zta, vpic, no_site], ["rta", "zta", "vpic", "none"]):
        
    h3k4me3_coverage = extract_coverage(H3K4ME3_str1_paths, regions=reg)
    h3k4me3_reac_coverage = extract_coverage(H3K4ME3_REAC_str1_paths, regions=reg)

    fig = plt.figure()
    ctl_normed = np.sum(h3k4me3_coverage, 0) / len(reg)
    induced_normed = np.sum(h3k4me3_reac_coverage, 0) / len(reg)
    plt.plot(range(3000,7000), ctl_normed[3000:7000])
    plt.fill_between(range(3000,7000), ctl_normed[3000:7000], alpha=.3)
    plt.plot(range(3000,7000), induced_normed[3000:7000])
    plt.fill_between(range(3000,7000), induced_normed[3000:7000], alpha=.3)    
    plt.xticks([5000])
    plt.xlim([3000,7000])
    plt.ylim([0.95*np.min([np.min(ctl_normed[3000:7000]), np.min(induced_normed[3000:7000])]), .03])
    plt.yticks([])
    plt.axvline([5000],c='k')
    plt.savefig(f'akata_summation_curve_h3k4me3{length}.{reg_name}.svg')
    plt.close('all')






for reg, reg_name in zip([rta, zta, vpic, no_site], ["rta", "zta", "vpic", "none"]):
            
    # h3k4me3_coverage = extract_coverage(H3K4ME3_str1_paths, regions=reg)
    h3k27me3_coverage = extract_coverage(H3K27ME3_str1_paths, regions=reg)
    h3k27ac_coverage = extract_coverage(H3K27AC_str1_paths, regions=reg)
    #h3k4me3_reac_coverage = extract_coverage(H3K4ME3_REAC_str1_paths, regions=reg)
    height = 8 * (len(reg) / len(regions))
    all_matrices = [h3k27me3_coverage, h3k27ac_coverage]
    for i,j in zip(all_matrices, ["h3k27me3", 'h3k27ac']):
       
        fig = plt.figure(figsize=(3, height))
        plt.imshow(i, cmap='Purples',aspect='auto', interpolation="gaussian", vmax=20)
        plt.xticks([length//2])
        plt.yticks([])
        plt.savefig(j + f'_denovo_histonecov_{length}.{reg_name}.vmax120.svg')
        plt.close('all')




canonical_bed_path = '/Users/nate/dnovo/beds/CANONICAL_TSS_from_mutu_ctl.noblacklist.bed'
canonical = import_bed(canonical_bed_path,col_sort=6)


h3k4me3_coverage = extract_coverage(H3K4ME3_str1_paths, regions=canonical)
h3k4me3_reac_coverage = extract_coverage(H3K4ME3_REAC_str1_paths, regions=canonical)
for i,j in zip([h3k4me3_coverage,h3k4me3_reac_coverage], ["h3k4me3", 'h3k4me3_reac']):
    
    fig = plt.figure(figsize=(3, 8))
    plt.imshow(i, cmap='Purples',aspect='auto', interpolation="gaussian", vmax=100)
    plt.xticks([length//2])
    plt.yticks([])
    plt.savefig(j + f'_denovo_histonecov_{length}.canonical.vmax120.svg')
    plt.close('all')


zta_chip_sites = '/Users/nate/dnovo/beds/Raji_Zta_ChIP_induced_summits_slop5000.no_blacklist.bed'
zta_chip = import_bed(zta_chip_sites, col_sort=4)

h3k4me3_coverage = extract_coverage(H3K4ME3_str1_paths, regions=zta_chip)
h3k4me3_reac_coverage = extract_coverage(H3K4ME3_REAC_str1_paths, regions=zta_chip)
for i,j in zip([h3k4me3_coverage,h3k4me3_reac_coverage], ["h3k4me3", 'h3k4me3_reac']):
    
    fig = plt.figure(figsize=(3, 8))
    plt.imshow(i, cmap='Purples',aspect='auto', interpolation="gaussian", vmax=100)
    plt.xticks([length//2])
    plt.yticks([])
    plt.savefig(j + f'_denovo_histonecov_{length}.zta_chip.vmax120.svg')
    plt.close('all')




h3k27me3_coverage = extract_coverage(H3K27ME3_str1_paths, regions=reg)
h3k27ac_coverage = extract_coverage(H3K27AC_str1_paths, regions=reg)


regions = []
with open(canonical_bed_path) as bed_handle:
    for line in bed_handle:
        regions.append(line.strip('\n').split('\t'))

random.shuffle(regions)
regions.sort(key=lambda x:float(x[5]))

        
h3k4me3_coverage = extract_coverage(H3K4ME3_str1_paths, regions=regions)
h3k27me3_coverage = extract_coverage(H3K27ME3_str1_paths, regions=regions)
h3k27ac_coverage = extract_coverage(H3K27AC_str1_paths, regions=regions)

    
all_matrices = [h3k4me3_coverage, h3k27me3_coverage, h3k27ac_coverage]
for i,j in zip(all_matrices, ["h3k4me3", "h3k27me3", "h3k27ac"]):
    fig = plt.figure(figsize=(3, 8))
    plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian", vmax=30)
    plt.xticks([length//2])
    plt.yticks([])
    plt.savefig(j + f'_canonical_histone_{length}.svg')
    plt.close('all')

