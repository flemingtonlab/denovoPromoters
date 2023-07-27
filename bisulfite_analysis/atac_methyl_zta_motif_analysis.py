import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

min_cov = 10

methyl_zta_site_path = '/Volumes/de_novo/ATAC/Mutu_12h_ATAC/bams/BINDetect_output_MC12_new_methyl_plus_all/meZta1_EBVmeZRE1/meZta1_EBVmeZRE1_overview.txt'
methyl_zta_atac = pd.read_table(methyl_zta_site_path)

ctls = ['C1', 'C2', 'C3', 'C4']
igs = ['Ig1', 'Ig2', 'Ig3', 'Ig4']

def get_path(prefix, chromosome):
    dirname = '/Volumes/de_novo/Bisulfite/Akata/coverage'
    return f'{dirname}/{prefix}.{prefix}_1_bismark_bt2_pe.deduplicated.bismark.cov_{chromosome}'


for sample in ctls:
    
    d = {}
    for chromosome in set(methyl_zta_atac['TFBS_chr']):
        filt_df = methyl_zta_atac[methyl_zta_atac['TFBS_chr']== chromosome]
        path = get_path(sample, chromosome)
        cov = pd.read_table(path, header=None)
        cov['total'] = cov[4] + cov[5]
        cov = cov[cov['total'] < min_cov]
        for row_id in filt_df.index:
            start, stop = filt_df.loc[row_id, ['TFBS_start', 'TFBS_end']]
            sites = cov[(cov[0]==chromosome) & (cov[1] > start) & (cov[1] < stop)]
            methylation = np.mean(sites[3])
            d[row_id] = methylation
    
    methyl_zta_atac[sample] = [d[i] for i in methyl_zta_atac.index]