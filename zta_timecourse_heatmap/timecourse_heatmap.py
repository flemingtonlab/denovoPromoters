import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from pathlib import Path
import re


direc = '/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure2/scripts/denovo_promoter_coverage_scripts/normalized_dn_coverage_timecourse_data/'
df6 = pd.read_table(direc + '6htime_summary.tsv.de.out', index_col=0)
df12 = pd.read_table(direc + '12htime_summary.tsv.de.out', index_col=0)
df24 = pd.read_table(direc + '24htime_summary.tsv.de.out', index_col=0)

df6 = df6[df6['padj']<.05]
df12 = df12[df12['padj']<.05]
df24 = df24[df24['padj'] <.05]

df12_only = df12.loc[set(df12.index) - set(df6.index)]
df24_only = df24.loc[set(df24.index) - set(df12_only.index)]

df6 = df6.sample(frac=1)
df12_only = df12_only.sample(frac=1)
df24_only = df24_only.sample(frac=1)

direc = '/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure2/scripts/denovo_promoter_coverage_scripts/normalized_dn_coverage_timecourse_data/'
n = pd.read_table(direc + 'concat_normalized_dn_expression300bp_log_fc_over_meanctl.tsv',sep='\t',index_col=0)
l = []
for i in [df6, df12_only, df24_only]:
    for key in i.index:
        if key not in l:
            l.append(key)
n2 = n.loc[l]


working_dir = "/Users/nate/Documents/Projects/De_Novo_Promoter/Manuscript/Figures/Figure2/starting_files/"

fa40_25 = {}
with open(working_dir + 'de_novo_tss_40to25bp_upstream.fa') as infile:
    for line in infile:
        line = line.strip('\n')
        if line[0] == '>':
            keya = line.replace('>','').split('(')[0]
        else:
            fa40_25[keya]= line
fa200_0 = {}
with open(working_dir + 'de_novo_tss_200bp_upstream.fa') as infile:
    for line in infile:
        line = line.strip('\n')
        if line[0] == '>':
            keya = line.replace('>','').split('(')[0]
        else:
            fa200_0[keya]= line

df = pd.read_table('/Users/nate/Downloads/motifs_in_de_novo_promoters/4_vPIC_plus_Zta/4_TATTAAA_TATTTAA_and_Zta_ChIP_seq_Hammerschmidt/TATTAAA_plus_TATTTAA_and_Zta_ChIP_Hammerschmidt_output/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed_TATTAAA_plus_TATTTAA.bed.binary_motif_info.bed_Raji_Zta_ChIP_induced_pooledReps_summits.bed.no_EBV_plus_strand.bed.binary_motif_info.tsv', index_col=0)
df = df.set_index('DN prom name')


n2['promoter_40_25'] = n2.index.map(lambda h:fa40_25[h])
n2['promoter_200_0'] = n2.index.map(lambda h:fa200_0[h])
bcrf1_motif = "TATT[TA]AA"
rta_motif_fw = "G[ACTG]CC[ACGT]{8,10}GG[ACGT]G" 
rta_motif_rev= "C[ACGT]CC[ACGT]{8,10}GG[ACGT]C"
rta_prog = re.compile(f'{rta_motif_fw}|{rta_motif_rev}')

bcrf1_prog = re.compile(bcrf1_motif)
n2['bcrf1_motif'] = n2.promoter_40_25.map(lambda y:bcrf1_prog.search(y) is not None)
n2['rta_motif'] = n2.promoter_200_0.map(lambda y:rta_prog.search(y) is not None)
# zta = pd.read_table(working_dir + 'de_novo_promoters_all_Mutu_Zta.corrected.slopped200l25r.ZTA_ChIP_intersect.bed', index_col=3, header=None)
# zta = zta.reset_index().drop_duplicates(3).set_index(3)
# n2['zta'] = zta[4]

#n2['bcrf1_motif'] = df['TATTAAA_plus_TATTTAA']
n2['zta'] = df['Zta_ChIP_Hammerschmidt']

n2_24h = n2.loc[df24_only.index]
n2_12h = n2.loc[df12_only.index]
n2_6h = n2.loc[df6.index]

zta_24 = ['blue' if i >0  else 'white' for i in n2_24h['zta']]
rta_24 = ['purple' if i is True else 'white' for i in n2_24h['rta_motif']]
bcrf1_24 = ['green' if i >0 else 'white' for i in n2_24h['bcrf1_motif']]

zta_12 = ['blue' if i >0  else 'white' for i in n2_12h['zta']]
rta_12 = ['purple' if i is True else 'white' for i in n2_12h['rta_motif']]
bcrf1_12 = ['green' if i >0 else 'white' for i in n2_12h['bcrf1_motif']]

zta_6 = ['blue' if i >0  else 'white' for i in n2_6h['zta']]
rta_6 = ['purple' if i is True else 'white' for i in n2_6h['rta_motif']]
bcrf1_6 = ['green' if i >0 else 'white' for i in n2_6h['bcrf1_motif']]


len24h = len(n2_24h.index)
len12h = len(n2_12h.index)
len6h = len(n2_6h.index)

max_height = np.max([len24h, len12h, len6h])
heightfrac24 = len24h / max_height
heightfrac12 = len12h / max_height
heightfrac6 = len6h / max_height


cg = sns.clustermap(n2_6h.iloc[:,:12],figsize=(20,20*heightfrac6), row_cluster=False, col_cluster=False, cmap='Reds',vmin=0,vmax=13, row_colors=[ rta_6,zta_6,bcrf1_6,],yticklabels=False, xticklabels=False)
cg.ax_row_dendrogram.set_visible(False)
cg.cax.set_visible(False)

plt.savefig("timecourse_heatmap_6h.svg")
plt.savefig("timecourse_heatmap_6h.png")

plt.close('all')

cg = sns.clustermap(n2_12h.iloc[:,:12],figsize=(20,20*heightfrac12), row_cluster=False, col_cluster=False, cmap='Reds',vmin=0,vmax=13, row_colors=[rta_12,zta_12,bcrf1_12,],yticklabels=False, xticklabels=False)
cg.ax_row_dendrogram.set_visible(False)
cg.cax.set_visible(False)

plt.savefig("timecourse_heatmap_12h.svg")
plt.savefig("timecourse_heatmap_12h.png")

plt.close('all')

cg = sns.clustermap(n2_24h.iloc[:,:12],figsize=(20,20*heightfrac24), row_cluster=False, col_cluster=False, cmap='Reds',vmin=0,vmax=13, row_colors=[rta_24,zta_24,bcrf1_24,],yticklabels=False, xticklabels=False)
cg.ax_row_dendrogram.set_visible(False)
cg.cax.set_visible(False)
plt.savefig("timecourse_heatmap_24h.svg")
plt.savefig("timecourse_heatmap_24h.png")

plt.close('all')
