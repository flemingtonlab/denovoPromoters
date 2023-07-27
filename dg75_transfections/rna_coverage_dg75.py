#!/usr/bin/env python

import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans

'''
6/29/23
1. Erik's Perl script was used to calculate denovo transcript abundance in DG75 cells
2. Values (rta vs ctl or zta vs ctl (n=4)) were input into deseq2
3. Sites with p<.05 for either zta or rta were analyzed using this script

'''
bed_path = '/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure2/starting_files/sites_induced_by_either_rta_or_zta.bed'

akata_uninduced_str1_paths = [
    '/Users/nate/dnovo/rna/akata_bcr/24hrUninduced_1.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrUninduced_2.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrUninduced_3.Aligned.out.sorted.str1.bw',
    ]

akata_uninduced_str2_paths = [
    '/Users/nate/dnovo/rna/akata_bcr/24hrUninduced_1.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrUninduced_2.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrUninduced_3.Aligned.out.sorted.str2.bw'
    ]
# human + ebv
# akata_uninduced_total_coverages =[64077140171, 43262093648, 54603139641]
akata_uninduced_total_coverages = [5.090669e+10, 3.769807e+10, 4.734046e+10]
akata_induced_str1_paths = [
    '/Users/nate/dnovo/rna/akata_bcr/24hrInduced_1.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrInduced_2.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrInduced_3.Aligned.out.sorted.str1.bw'
    ]

akata_induced_str2_paths = [
    '/Users/nate/dnovo/rna/akata_bcr/24hrInduced_1.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrInduced_2.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/akata_bcr/24hrInduced_3.Aligned.out.sorted.str2.bw'
    ]

# human + ebv
# akata_induced_total_coverages =[35194259069, 27213959966, 27111502443]
akata_induced_total_coverages = [1.108861e+10, 9.385569e+09, 9.250357e+09]

mutu_ctl_str1_paths = [
    '/Users/nate/dnovo/rna/mutu_zta/MC1.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MC2.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MC3.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MC4.Aligned.out.sorted.str1.bw'
]

mutu_ctl_str2_paths = [
    '/Users/nate/dnovo/rna/mutu_zta/MC1.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MC2.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MC3.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MC4.Aligned.out.sorted.str2.bw'
]

#human+ebv
#mutu_ctl_total_coverages =[37260586565, 33417299733, 36167651053, 35774470154]
mutu_ctl_total_coverages = [3.492805e+10, 3.139268e+10, 3.420245e+10, 3.390367e+10]

mutu_zta_str1_paths = [
    '/Users/nate/dnovo/rna/mutu_zta/MZ1.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MZ3.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MZ4.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MZ5.Aligned.out.sorted.str1.bw',

]

mutu_zta_str2_paths = [
    '/Users/nate/dnovo/rna/mutu_zta/MZ1.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MZ3.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MZ4.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_zta/MZ5.Aligned.out.sorted.str2.bw',

]
#human + ebv
#mutu_zta_total_coverages =[30409351621, 39773538572, 33328590766, 31399955669]
mutu_zta_total_coverages = [1.545093e+10, 2.021701e+10, 1.703477e+10, 1.518691e+10]


mutu_uninduced_str1_paths = [
    '/Users/nate/dnovo/rna/mutu_bcr/MU1.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MU2.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MU3.Aligned.out.sorted.str1.bw',
]

mutu_uninduced_str2_paths = [
    '/Users/nate/dnovo/rna/mutu_bcr/MU1.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MU2.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MU3.Aligned.out.sorted.str2.bw',
]

#human + ebv
#mutu_uninduced_total_coverages =[27521288759, 27986327116, 27346255520]
mutu_uninduced_total_coverages = [2.608932e+10, 2.661923e+10, 2.612005e+10]

mutu_induced_str1_paths = [
    '/Users/nate/dnovo/rna/mutu_bcr/MI1.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MI2.Aligned.out.sorted.str1.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MI3.Aligned.out.sorted.str1.bw',
]

mutu_induced_str2_paths = [
    '/Users/nate/dnovo/rna/mutu_bcr/MI1.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MI2.Aligned.out.sorted.str2.bw',
    '/Users/nate/dnovo/rna/mutu_bcr/MI3.Aligned.out.sorted.str2.bw',
]

# human + ebv
# mutu_induced_total_coverages =[30948732341, 30923000690, 24373351477]
mutu_induced_total_coverages = [1.444843e+10, 1.423205e+10, 1.165600e+10]

dg75_ctl_str1_paths = [
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C1.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C2.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C3.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C4.Aligned.out.sorted.str1.bw'
]
dg75_ctl_str2_paths = [
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C1.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C2.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C3.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_C4.Aligned.out.sorted.str2.bw'
]
dg75_zta_str1_paths = [
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z1.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z2.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z3.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z4.Aligned.out.sorted.str1.bw'
]
dg75_zta_str2_paths = [
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z1.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z2.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z3.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_Z4.Aligned.out.sorted.str2.bw'
]


dg75_rta_str1_paths = [
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R1.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R2.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R3.Aligned.out.sorted.str1.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R4.Aligned.out.sorted.str1.bw'
]
dg75_rta_str2_paths = [
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R1.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R2.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R3.Aligned.out.sorted.str2.bw',
    '/Volumes/de_novo/RNA/DG75_Zta_Rta/bigwigs/DG75_R4.Aligned.out.sorted.str2.bw'
]


stranded = True

upstream_bases = 1000
downstream_bases = 1000
length = upstream_bases + downstream_bases

regions = []
with open(bed_path) as bed_handle:
    for line in bed_handle:
        regions.append(line.strip('\n').split('\t'))

def extract_coverage(str1_paths, regions):
    
    m = np.zeros([len(regions), length])

    for str1_path in str1_paths:
        str1 = pyBigWig.open(str1_path)
        str2 = pyBigWig.open(str1_path.replace('str1','str2'))


        for ind, region in enumerate(regions):


            start = int(float(region[1]))
            stop = int(float(region[2]))
            middle = (start + stop) //2
            start = middle - upstream_bases
            stop = middle + downstream_bases
            strand = region[5]

            if strand == '+':
                vals = str2.values(region[0], start, stop)
            else:
                vals = str1.values(region[0], start, stop)
                vals = np.flip(vals)
            vals = np.nan_to_num(np.array(vals), nan=0)
            vals = np.abs(vals)
            #vals = vals/np.max(vals+.001)

            m[ind] +=  vals
        print(str1_path, "done")

    m/=len(str1_paths)
    return m 


def get_row_cov_max(matrices_list, min_divisor=5):
    rows = matrices_list[0].shape[0]
    number_of_matrices = len(matrices_list)
    maxes = np.zeros([rows, number_of_matrices + 1])
    for col, matr in enumerate(matrices_list):
        maxes[:, col] = np.max(matr, 1)
    maxes[number_of_matrices] = min_divisor
    return np.max(maxes, 1)


akata_uninduced_coverage = extract_coverage(akata_uninduced_str1_paths, regions=regions)
akata_induced_coverage = extract_coverage(akata_induced_str1_paths, regions=regions)
mutu_uninduced_coverage = extract_coverage(mutu_uninduced_str1_paths, regions=regions)
mutu_induced_coverage = extract_coverage(mutu_induced_str1_paths, regions=regions)
mutu_ctl_coverage = extract_coverage(mutu_ctl_str1_paths, regions=regions)
mutu_zta_coverage = extract_coverage(mutu_zta_str1_paths,  regions=regions)
dg75_ctl_coverage = extract_coverage(dg75_ctl_str1_paths,  regions=regions)
dg75_zta_coverage = extract_coverage(dg75_zta_str1_paths, regions=regions)
dg75_rta_coverage = extract_coverage(dg75_rta_str1_paths,  regions=regions)


all_matrices = [
    # akata_uninduced_coverage,
    # akata_induced_coverage,
    # mutu_uninduced_coverage,
    # mutu_induced_coverage,
    # mutu_ctl_coverage,
    # mutu_zta_coverage,
    dg75_ctl_coverage,
    dg75_zta_coverage,
    dg75_rta_coverage
    ]


max_all = get_row_cov_max(all_matrices, .5)
for ind in range(len(all_matrices)):
    all_matrices[ind] /= max_all[:, None]



for i,j in zip(all_matrices, ['dg_c','dg_z','dg_r']):#["ak_un", "ak_in", "mu_un", "mu_in", "mu_c","mu_z",'dg_c','dg_z','dg_r']):
    fig = plt.figure(figsize=(3,2))
    plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian",vmax=.6)
    plt.xticks([1000])
    plt.yticks([])
    plt.savefig(j + 'dg75_sites.1k_bsort.svg')
    plt.close('all')

all_matrices = [
    mutu_ctl_coverage,
    mutu_zta_coverage
]
max_all = get_row_cov_max(all_matrices, .5)
for ind in range(len(all_matrices)):
    all_matrices[ind] /= max_all[:, None]

for i,j in zip(all_matrices, ['mu_c','mu_z']):#["ak_un", "ak_in", "mu_un", "mu_in", "mu_c","mu_z",'dg_c','dg_z','dg_r']):
    fig = plt.figure(figsize=(3,2))
    plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian",vmax=.6)
    plt.xticks([1000])
    plt.yticks([])
    plt.savefig(j + 'dg75_sites.1k_bsort.svg')
    plt.close('all')


# Did not use kmeans for figure
kmeans = KMeans(n_clusters=10).fit(all_matrices[1])
a = pd.DataFrame(kmeans.labels_)

all_matrices_km = []
for i in all_matrices:
    all_matrices_km.append(i[np.array(a.sort_values(0).index)])

for i,j in zip(all_matrices_km, ["ak_un", "ak_in", "mu_un", "mu_in", "mu_c","mu_z"]):
    fig = plt.figure(figsize=(3,8))
    plt.imshow(i, cmap='Reds',aspect='auto', interpolation="gaussian",vmax=.4)
    plt.xticks([1000])
    plt.yticks([])
    plt.savefig(j + '1k.kmeans_b.svg')
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
plot_summution_curves(all_matrices, [ "mu_c","mu_z"])