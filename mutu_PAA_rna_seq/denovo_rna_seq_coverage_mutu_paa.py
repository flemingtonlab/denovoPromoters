#!/usr/bin/env python

import sys
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import re
import glob

# PAA experiment done in Mutu cells
bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_all_Mutu_Zta.bed'

CnoP_str1 = glob.glob('/Volumes/de_novo/RNA/Mutu_PAA/bigwigs/MC_*.Aligned.out.sorted.str1.bw')
CnoP_str2 = [i.replace('str1','str2') for i in CnoP_str1]

ZnoP_str1 = glob.glob('/Volumes/de_novo/RNA/Mutu_PAA/bigwigs/MZ_[1-9].Aligned.out.sorted.str1.bw')
ZnoP_str2 = [i.replace('str1','str2') for i in ZnoP_str1]

ZP_str1 = glob.glob('/Volumes/de_novo/RNA/Mutu_PAA/bigwigs/MZ_P_*.Aligned.out.sorted.str1.bw')
ZP_str2 = [i.replace('str1','str2') for i in ZP_str1]

# total_reads=$(samtools idxstats $bam |grep -v "EBV" | cut -f3 |paste -s -d+ - |bc)
CnoP_coverages = [247340825, 223181897, 244408212, 248942510]
ZnoP_coverages = [147064659, 143462143, 140344108, 145032402]
ZP_coverages = [181457180, 193401417, 182974955, 179455703]

stranded = True
upstream_bases = 1000
downstream_bases = 1000
length = upstream_bases + downstream_bases

regions = []
with open(bed_path) as bed_handle:
    for line in bed_handle:
        regions.append(line.strip('\n').split('\t'))
regions.sort(key=lambda x:float(x[4]))


def extract_coverage(strand1_paths, strand2_paths, total_coverages, regions):
    
    m = np.zeros([len(regions), length])
    for str1_path, str2_path, total_coverage in zip(strand1_paths, strand2_paths, total_coverages):
        str1 = pyBigWig.open(str1_path)
        str2 = pyBigWig.open(str2_path)
        for ind, region in enumerate(regions):
            start = int(float(region[1]))
            stop = int(float(region[2])) 
            strand = region[5]
            middle = (start + stop) // 2 
            start = middle - upstream_bases
            stop = middle + downstream_bases
            if stranded and region[5] == '+':
                vals = str2.values(region[0], start, stop)
            elif stranded and region[5] == '-':
                vals = str1.values(region[0], start, stop)
            if strand == '-':
                vals = np.flip(vals)         
            vals = np.abs(vals)
            vals = np.nan_to_num(np.array(vals), nan=0.001)
            m[ind] += (1000000 * vals)/ total_coverage       
        print(str1_path, "done")
    return m / len(strand1_paths)


def get_row_cov_max(matrices_list, min_divisor=10):
    rows = matrices_list[0].shape[0]
    number_of_matrices = len(matrices_list)
    maxes = np.zeros([rows, number_of_matrices + 1])
    for col, matr in enumerate(matrices_list):
        maxes[:, col] = np.max(matr, 1)
    maxes[number_of_matrices] = min_divisor
    return np.max(maxes, 1)


def plot_heatmap(matrix, name, colormap='Reds',height=8):
    fig = plt.figure(figsize=(3, height))
    plt.imshow(matrix, cmap=colormap, aspect='auto', interpolation="gaussian", vmax=.6)
    plt.xticks([])
    plt.yticks([])
    plt.savefig(name + '.heatmap.svg')
    plt.close('all') 


def plot_summution_curves(matrices_list, matrix_names):
    max_val = 0
    for matr in matrices_list:
        matr_sum = np.sum(matr, 0)
        matr_max = np.max(matr_sum)
        if matr_max > max_val:
            max_val = matr_max

    for matr, matr_name in zip(matrices_list, matrix_names):
        fig = plt.figure()
        plt.plot(range(matr[0].shape[0]), np.sum(matr, 0), c='k') 
        plt.ylim([0, max_val])  
        plt.xlim([0, matr[0].shape[0]])
        plt.yticks([])
        plt.xticks([])
        plt.savefig(matr_name + 'summation_curve.svg') 
        plt.close('all')


zta_vpic_df = pd.read_table('/Users/nate/Downloads/de_novo_promoter_paper/motifs_in_de_novo_promoters/4_vPIC_plus_Zta/4_TATTAAA_TATTTAA_and_Zta_ChIP_seq_Hammerschmidt/TATTAAA_plus_TATTTAA_and_Zta_ChIP_Hammerschmidt_output/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed_TATTAAA_plus_TATTTAA.bed.binary_motif_info.bed_Raji_Zta_ChIP_induced_pooledReps_summits.bed.no_EBV_plus_strand.bed.binary_motif_info.tsv', index_col=0)
zta_vpic_df = zta_vpic_df.set_index('DN prom name')
zta = zta_vpic_df[zta_vpic_df['Zta_ChIP_Hammerschmidt'] == 1]
zta = set(zta.index)
zta_regions = [i for i in regions if i[3] in zta]

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

fa200_0 = {}
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
rta_regions = [i for i in regions if rta_prog.search(fa200_0[i[3]]) is not None]
vpic_regions = [i for i in regions if bcrf1_prog.search(fa40_25[i[3]]) is not None]
all_three = rta_regions+zta_regions+vpic_regions
all_three = [i[3] for i in all_three]
l = []
two_plus_sites = []
for i in all_three:
    if i in l:
        two_plus_sites.append(i)
    l.append(i)
no_site = [i for i in regions if i[3] not in all_three]

vpic_regions = [i for i in vpic_regions if i not in two_plus_sites]
zta_regions = [i for i in zta_regions if i not in two_plus_sites]
rta_regions = [i for i in rta_regions if i not in two_plus_sites]

vpic_C_cov = extract_coverage(CnoP_str1, CnoP_str2, CnoP_coverages, vpic_regions)
vpic_Z_cov = extract_coverage(ZnoP_str1, ZnoP_str2, ZnoP_coverages, vpic_regions)
vpic_ZP_cov = extract_coverage(ZP_str1, ZP_str2, ZP_coverages, vpic_regions)
all_matrices = [vpic_C_cov, vpic_Z_cov, vpic_ZP_cov]
maxes = get_row_cov_max(all_matrices, 5)
C_cov_m = vpic_C_cov / maxes[:, None]
Z_cov_m = vpic_Z_cov / maxes[:, None]
ZP_cov_m = vpic_ZP_cov / maxes[:, None]
plot_heatmap(C_cov_m, "CC.vpic.rna.coverage", 'Reds', 8*(len(vpic_regions)/len(regions)))
plot_heatmap(Z_cov_m, "ZC.vpic.rna.coverage", 'Reds', 8*(len(vpic_regions)/len(regions)))
plot_heatmap(ZP_cov_m, "ZPaa.vpic.rna.coverage", 'Reds', 8*(len(vpic_regions)/len(regions)))
plot_summution_curves([C_cov_m, Z_cov_m, ZP_cov_m], ["C_cov_m_vpic", "Z_cov_m_vpic", "ZP_cov_m_vpic"])

zta_C_cov = extract_coverage(CnoP_str1, CnoP_str2, CnoP_coverages, zta_regions)
zta_Z_cov = extract_coverage(ZnoP_str1, ZnoP_str2, ZnoP_coverages, zta_regions)
zta_ZP_cov = extract_coverage(ZP_str1, ZP_str2, ZP_coverages, zta_regions)
all_matrices = [zta_C_cov, zta_Z_cov, zta_ZP_cov]
maxes = get_row_cov_max(all_matrices, 5)
C_cov_m = zta_C_cov / maxes[:, None]
Z_cov_m = zta_Z_cov / maxes[:, None]
ZP_cov_m = zta_ZP_cov / maxes[:, None]
plot_heatmap(C_cov_m, "CC.zta.rna.coverage", 'Reds', 8*(len(zta_regions)/len(regions)))
plot_heatmap(Z_cov_m, "ZC.zta.rna.coverage", 'Reds', 8*(len(zta_regions)/len(regions)))
plot_heatmap(ZP_cov_m, "ZPaa.zta.rna.coverage", 'Reds', 8*(len(zta_regions)/len(regions)))
plot_summution_curves([C_cov_m, Z_cov_m, ZP_cov_m], ["C_cov_m_zta", "Z_cov_m_zta", "ZP_cov_m_zta"])

rta_C_cov = extract_coverage(CnoP_str1, CnoP_str2, CnoP_coverages, rta_regions)
rta_Z_cov = extract_coverage(ZnoP_str1, ZnoP_str2, ZnoP_coverages, rta_regions)
rta_ZP_cov = extract_coverage(ZP_str1, ZP_str2, ZP_coverages, rta_regions)
all_matrices = [rta_C_cov, rta_Z_cov, rta_ZP_cov]
maxes = get_row_cov_max(all_matrices, 5)
C_cov_m = rta_C_cov / maxes[:, None]
Z_cov_m = rta_Z_cov / maxes[:, None]
ZP_cov_m = rta_ZP_cov / maxes[:, None]
plot_heatmap(C_cov_m, "CC.rta.rna.coverage", 'Reds', 8*(len(rta_regions)/len(regions)))
plot_heatmap(Z_cov_m, "ZC.rta.rna.coverage", 'Reds', 8*(len(rta_regions)/len(regions)))
plot_heatmap(ZP_cov_m, "ZPaa.rta.rna.coverage", 'Reds', 8*(len(rta_regions)/len(regions)))
plot_summution_curves([C_cov_m, Z_cov_m, ZP_cov_m], ["C_cov_m_rta", "Z_cov_m_rta", "ZP_cov_m_rta"])

no_site_C_cov = extract_coverage(CnoP_str1, CnoP_str2, CnoP_coverages, no_site)
no_site_Z_cov = extract_coverage(ZnoP_str1, ZnoP_str2, ZnoP_coverages, no_site)
no_site_ZP_cov = extract_coverage(ZP_str1, ZP_str2, ZP_coverages, no_site)
all_matrices = [no_site_C_cov, no_site_Z_cov, no_site_ZP_cov]
maxes = get_row_cov_max(all_matrices, 5)
C_cov_m = no_site_C_cov / maxes[:, None]
Z_cov_m = no_site_Z_cov / maxes[:, None]
ZP_cov_m = no_site_ZP_cov / maxes[:, None]
plot_heatmap(C_cov_m, "CC.no_site.rna.coverage", 'Reds', 8*(len(no_site)/len(regions)))
plot_heatmap(Z_cov_m, "ZC.no_site.rna.coverage", 'Reds', 8*(len(no_site)/len(regions)))
plot_heatmap(ZP_cov_m, "ZPaa.no_site.rna.coverage", 'Reds', 8*(len(no_site)/len(regions)))
plot_summution_curves([C_cov_m, Z_cov_m, ZP_cov_m], ["C_cov_m_nosite", "Z_cov_m_nosite", "ZP_cov_m_nosite"])