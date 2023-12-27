
import sys
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import random
import re
import glob

# Plot created July 7, 2023, modified Dec 22, 2023
# Any de novo site that had overlap with a blacklisted region (encode) was removed (using the -5000 to +5000 region)
bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.noblacklist.bed'

bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.noblacklist.bed'
canonical_bed_path = '/Users/nate/dnovo/beds/TSS_from_mutu_ctl2.bed'

akata_uninduced_str1_paths = ['/Volumes/de_novo/ATAC/Akata_ATAC/bams/AkataC24h.footprints.bw']
akata_induced_str1_paths = ['/Volumes/de_novo/ATAC/Akata_ATAC/bams/AkataIg24h.footprints.bw']

mutu_ctl_str1_paths = ['/Volumes/de_novo/ATAC/Mutu_ATAC/bams/MC24.footprints.bw']
mutu_zta_str1_paths = ['/Volumes/de_novo/ATAC/Mutu_ATAC/bams/MZ24.footprints.bw']

upstream_bases = 5000
downstream_bases = 5000
length = upstream_bases + downstream_bases

regions = []
with open(bed_path) as bed_handle:
    for line in bed_handle:
        regions.append(line.strip('\n').split('\t'))

random.shuffle(regions)
regions.sort(key=lambda x:float(x[4]))


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
vpic_regions = [i for i in vpic_regions if i[3] not in two_plus_sites]
zta_regions = [i for i in zta_regions if i[3] not in two_plus_sites]
rta_regions = [i for i in rta_regions if i[3] not in two_plus_sites]

for reg, reg_name in zip([rta_regions, zta_regions, vpic_regions, no_site], ["rta", "zta", "vpic", "none"]):
        
    akata_uninduced_coverage = extract_coverage(akata_uninduced_str1_paths, regions=reg)
    akata_induced_coverage = extract_coverage(akata_induced_str1_paths, regions=reg)
    height = 8 * (len(reg) / 12000)
    all_matrices = [akata_uninduced_coverage, akata_induced_coverage]
    for i,j in zip(all_matrices, ["ak_un", "ak_in"]):
        fig = plt.figure(figsize=(3, height))
        plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian", vmax=.5)
        plt.colorbar()

        plt.xticks([length//2])
        plt.yticks([])
        plt.savefig(j + f'_denovo_noblacklist_atac_akata_{length}.{reg_name}.svg')
        plt.close('all')

for reg, reg_name in zip([rta_regions, zta_regions, vpic_regions, no_site], ["rta", "zta", "vpic", "none"]):
        
    mutu_ctl_coverage = extract_coverage(mutu_ctl_str1_paths, regions=reg)
    mutu_zta_coverage = extract_coverage(mutu_zta_str1_paths, regions=reg)

    height = 8 * (len(reg) / 12000)
    all_matrices = [mutu_ctl_coverage, mutu_zta_coverage]
    for i,j in zip(all_matrices, ["mut_c", "mut_z"]):
        fig = plt.figure(figsize=(3, height))
        plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian", vmax=.5)
        plt.xticks([length//2])
        plt.yticks([])
        plt.savefig(j + f'_denovo_noblacklist_atac_mutu_{length}.{reg_name}.svg')
        plt.close('all')

for reg, reg_name in zip([rta_regions, zta_regions, vpic_regions, no_site], ["rta", "zta", "vpic", "none"]):
        
    mutu_ctl_coverage = extract_coverage(mutu_ctl_str1_paths, regions=reg)
    mutu_zta_coverage = extract_coverage(mutu_zta_str1_paths, regions=reg)

    fig = plt.figure()
    muctl_normed = np.sum(mutu_ctl_coverage, 0) / len(reg)
    muzta_normed = np.sum(mutu_zta_coverage, 0) / len(reg)
    plt.plot(range(3000,7000), muctl_normed[3000:7000])
    plt.fill_between(range(3000,7000), muctl_normed[3000:7000], alpha=.3)
    plt.plot(range(3000,7000), muzta_normed[3000:7000])
    plt.fill_between(range(3000,7000), muzta_normed[3000:7000], alpha=.3)    
    plt.xticks([5000])
    plt.xlim([3000,7000])
    plt.ylim([0.95*np.min([np.min(muctl_normed[3000:7000]), np.min(muzta_normed[3000:7000])]), 0.4])
    plt.yticks([])
    plt.axvline([5000],c='k')
    plt.savefig(f'mutu_summation_curve_atac_noblacklist_{length}.{reg_name}.svg')
    plt.close('all')


for reg, reg_name in zip([rta_regions, zta_regions, vpic_regions, no_site], ["rta", "zta", "vpic", "none"]):
        
    akata_uninduced_coverage = extract_coverage(akata_uninduced_str1_paths, regions=reg)
    akata_induced_coverage = extract_coverage(akata_induced_str1_paths, regions=reg)
    fig = plt.figure()
    akun_normed = np.sum(akata_uninduced_coverage, 0) / len(reg)
    akin_normed = np.sum(akata_induced_coverage, 0) / len(reg)
    plt.plot(range(3000,7000), akun_normed[3000:7000])
    plt.fill_between(range(3000,7000), akun_normed[3000:7000], alpha=.3)
    plt.plot(range(3000,7000), akin_normed[3000:7000])
    plt.fill_between(range(3000,7000), akin_normed[3000:7000], alpha=.3)    
    plt.xticks([5000])
    plt.xlim([3000,7000])
    plt.ylim([0.95*np.min([np.min(akun_normed[3000:7000]), np.min(akin_normed[3000:7000])]), 0.40])
    plt.yticks([])
    plt.axvline([5000],c='k')
    plt.savefig(f'akata_summation_curve_atac_noblacklist_{length}.{reg_name}.svg')
    plt.close('all')
