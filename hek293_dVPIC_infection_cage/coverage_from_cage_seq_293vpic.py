#!/usr/bin/env python
import sys
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
import random
import re

stranded = True
upstream_bases = 2000
downstream_bases = 2000
length = upstream_bases + downstream_bases
denovo_bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed'


# HEK293 experiment, Djavadia et al. https://doi.org/10.1371/journal.ppat.1007114 
wt_induced_positive_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164144.fastq.gzSignal.Unique.str1.out.bw']
dOriLyt_induced_positive_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164145.fastq.gzSignal.Unique.str1.out.bw']
dBALF2_uninduced_positive_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164146.fastq.gzSignal.Unique.str1.out.bw']
dBALF2_induced_positive_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164147.fastq.gzSignal.Unique.str1.out.bw']
dBALF2_induced_transcomp_positive_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164148.fastq.gzSignal.Unique.str1.out.bw']
dBDLF4_induced_positive_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164149.fastq.gzSignal.Unique.str1.out.bw']
dBDLF4_induced_transcomp_positive_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164150.fastq.gzSignal.Unique.str1.out.bw']

wt_induced_negative_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164144.fastq.gzSignal.Unique.str2.out.bw']
dOriLyt_induced_negative_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164145.fastq.gzSignal.Unique.str2.out.bw']
dBALF2_uninduced_negative_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164146.fastq.gzSignal.Unique.str2.out.bw']
dBALF2_induced_negative_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164147.fastq.gzSignal.Unique.str2.out.bw']
dBALF2_induced_transcomp_negative_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164148.fastq.gzSignal.Unique.str2.out.bw']
dBDLF4_induced_negative_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164149.fastq.gzSignal.Unique.str2.out.bw']
dBDLF4_induced_transcomp_negative_strand = ['/Volumes/de_novo/uwmadison/CAGE/bigwig/SRR7164150.fastq.gzSignal.Unique.str2.out.bw']

# total_reads=$(samtools idxstats $bam |grep -v "EBV" | cut -f3 |paste -s -d+ - |bc)
wt_induced_total_reads = [36522781]
dOriLyt_induced_total_reads = [38272323]
dBALF2_uninduced_total_reads = [40722382]
dBALF2_induced_total_reads = [56376131]
dBALF2_induced_transcomp_total_reads = [45370505]
dBDLF4_induced_total_reads = [35358359]
dBDLF4_induced_transcomp_total_reads = [36913648]


def extract_coverage(strand1_paths, strand2_paths, total_coverages, regions, antisense=False):
    
    m = np.zeros([len(regions), length])
    for str1_path, str2_path, total_coverage in zip(strand1_paths, strand2_paths, total_coverages):
        if antisense:
            str1 = pyBigWig.open(str2_path)
            str2 = pyBigWig.open(str1_path)
        else:
            str1 = pyBigWig.open(str1_path)
            str2 = pyBigWig.open(str2_path)

        for ind, region in enumerate(regions):
            start = stop = int(float(region[1]))
            strand = region[5]
            if strand == '+':
                vals = str2.values(region[0], start-upstream_bases, stop+downstream_bases)
            else:
                vals = str1.values(region[0], start-downstream_bases, stop+upstream_bases)
                vals = np.flip(vals)                        
            vals = np.abs(vals)
            vals = np.nan_to_num(np.array(vals), nan=0)
            m[ind] += (1000000 * vals) / total_coverage
        print(str1_path, "done")
    return m / len(strand1_paths)


def plot_heatmap(matrix, name, colormap='Reds', height=8):
    fig = plt.figure(figsize=(3, height))
    plt.imshow(matrix, cmap=colormap, aspect='auto', interpolation="gaussian", vmax=.0002)
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


denovo = []
with open(denovo_bed_path) as bed_handle:
    for line in bed_handle:
        denovo.append(line.strip('\n').split('\t'))
random.shuffle(denovo)
denovo.sort(key=lambda x:float(x[4]))

zta_vpic_df = pd.read_table('/Users/nate/Downloads/de_novo_promoter_paper/motifs_in_de_novo_promoters/4_vPIC_plus_Zta/4_TATTAAA_TATTTAA_and_Zta_ChIP_seq_Hammerschmidt/TATTAAA_plus_TATTTAA_and_Zta_ChIP_Hammerschmidt_output/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed_TATTAAA_plus_TATTTAA.bed.binary_motif_info.bed_Raji_Zta_ChIP_induced_pooledReps_summits.bed.no_EBV_plus_strand.bed.binary_motif_info.tsv', index_col=0)
zta_vpic_df = zta_vpic_df.set_index('DN prom name')
zta = zta_vpic_df[zta_vpic_df['Zta_ChIP_Hammerschmidt'] == 1]
zta = set(zta.index)
zta_regions = [i for i in denovo if i[3] in zta]


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
rta_regions = [i for i in denovo if rta_prog.search(fa200_0[i[3]]) is not None]

vpic_regions = [i for i in denovo if bcrf1_prog.search(fa40_25[i[3]]) is not None]
all_three = rta_regions + zta_regions + vpic_regions
all_three = [i[3] for i in all_three]

l = []
two_plus_sites = []
for i in all_three:
    if i in l:
        two_plus_sites.append(i)
    l.append(i)

no_site = [i for i in denovo if i[3] not in all_three]
vpic_regions = [i for i in vpic_regions if i not in two_plus_sites]
zta_regions = [i for i in zta_regions if i not in two_plus_sites]
rta_regions = [i for i in rta_regions if i not in two_plus_sites]


vpic_wt_induced_coverage = extract_coverage(wt_induced_negative_strand, wt_induced_positive_strand, wt_induced_total_reads,regions=vpic_regions)
vpic_dOriLyt_induced_coverage = extract_coverage(dOriLyt_induced_negative_strand, dOriLyt_induced_positive_strand, dOriLyt_induced_total_reads,regions=vpic_regions)
vpic_dBALF2_uninduced_coverage = extract_coverage(dBALF2_uninduced_negative_strand, dBALF2_uninduced_positive_strand, dBALF2_uninduced_total_reads,regions=vpic_regions)
vpic_dBALF2_induced_coverage = extract_coverage(dBALF2_induced_negative_strand, dBALF2_induced_positive_strand, dBALF2_induced_total_reads,regions=vpic_regions)
vpic_dBALF2_induced_transcomp_coverage = extract_coverage(dBALF2_induced_transcomp_negative_strand, dBALF2_induced_transcomp_positive_strand, dBALF2_induced_transcomp_total_reads,regions=vpic_regions)
vpic_dBDLF4_induced_coverage = extract_coverage(dBDLF4_induced_negative_strand, dBDLF4_induced_positive_strand, dBDLF4_induced_total_reads,regions=vpic_regions)
vpic_dBDLF4_induced_transcomp_coverage = extract_coverage(dBDLF4_induced_transcomp_negative_strand, dBDLF4_induced_transcomp_positive_strand, dBDLF4_induced_transcomp_total_reads,regions=vpic_regions)

plot_heatmap(vpic_wt_induced_coverage, "vpic_wt_induced_sense.denovo.CAGE.coverage", 'Reds', 8* (len(vpic_regions)/len(denovo)))
plot_heatmap(vpic_dOriLyt_induced_coverage, "vpic_dOriLyt_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(vpic_regions)/len(denovo)))
plot_heatmap(vpic_dBALF2_uninduced_coverage, "vpic_dBALF2_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(vpic_regions)/len(denovo)))
plot_heatmap(vpic_dBALF2_induced_coverage, "vpic_dBALF2_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(vpic_regions)/len(denovo)))
plot_heatmap(vpic_dBALF2_induced_transcomp_coverage, "vpic_dBALF2_induced_transcomp_sense.denovo.CAGE.coverage", 'Reds', 8*(len(vpic_regions)/len(denovo)))
plot_heatmap(vpic_dBDLF4_induced_coverage, "vpic_dBDLF4_induced_sense.denovo.CAGE.coverage", "Reds", 8*(len(vpic_regions)/len(denovo)))
plot_heatmap(vpic_dBDLF4_induced_transcomp_coverage, "vpic_dBDLF4_induced_transcomp_sense.denovo.CAGE.coverage", "Reds", 8*(len(vpic_regions)/len(denovo)))


zta_wt_induced_coverage = extract_coverage(wt_induced_negative_strand, wt_induced_positive_strand, wt_induced_total_reads,regions=zta_regions)
zta_dOriLyt_induced_coverage = extract_coverage(dOriLyt_induced_negative_strand, dOriLyt_induced_positive_strand, dOriLyt_induced_total_reads,regions=zta_regions)
zta_dBALF2_uninduced_coverage = extract_coverage(dBALF2_uninduced_negative_strand, dBALF2_uninduced_positive_strand, dBALF2_uninduced_total_reads,regions=zta_regions)
zta_dBALF2_induced_coverage = extract_coverage(dBALF2_induced_negative_strand, dBALF2_induced_positive_strand, dBALF2_induced_total_reads,regions=zta_regions)
zta_dBALF2_induced_transcomp_coverage = extract_coverage(dBALF2_induced_transcomp_negative_strand, dBALF2_induced_transcomp_positive_strand, dBALF2_induced_transcomp_total_reads,regions=zta_regions)
zta_dBDLF4_induced_coverage = extract_coverage(dBDLF4_induced_negative_strand, dBDLF4_induced_positive_strand, dBDLF4_induced_total_reads,regions=zta_regions)
zta_dBDLF4_induced_transcomp_coverage = extract_coverage(dBDLF4_induced_transcomp_negative_strand, dBDLF4_induced_transcomp_positive_strand, dBDLF4_induced_transcomp_total_reads,regions=zta_regions)

plot_heatmap(zta_wt_induced_coverage, "zta_wt_induced_sense.denovo.CAGE.coverage", 'Reds', 8* (len(zta_regions)/len(denovo)))
plot_heatmap(zta_dOriLyt_induced_coverage, "zta_dOriLyt_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(zta_regions)/len(denovo)))
plot_heatmap(zta_dBALF2_uninduced_coverage, "zta_dBALF2_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(zta_regions)/len(denovo)))
plot_heatmap(zta_dBALF2_induced_coverage, "zta_dBALF2_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(zta_regions)/len(denovo)))
plot_heatmap(zta_dBALF2_induced_transcomp_coverage, "zta_dBALF2_induced_transcomp_sense.denovo.CAGE.coverage", 'Reds', 8*(len(zta_regions)/len(denovo)))
plot_heatmap(zta_dBDLF4_induced_coverage, "zta_dBDLF4_induced_sense.denovo.CAGE.coverage", "Reds", 8*(len(zta_regions)/len(denovo)))
plot_heatmap(zta_dBDLF4_induced_transcomp_coverage, "zta_dBDLF4_induced_transcomp_sense.denovo.CAGE.coverage", "Reds", 8*(len(zta_regions)/len(denovo)))


rta_wt_induced_coverage = extract_coverage(wt_induced_negative_strand, wt_induced_positive_strand, wt_induced_total_reads,regions=rta_regions)
rta_dOriLyt_induced_coverage = extract_coverage(dOriLyt_induced_negative_strand, dOriLyt_induced_positive_strand, dOriLyt_induced_total_reads,regions=rta_regions)
rta_dBALF2_uninduced_coverage = extract_coverage(dBALF2_uninduced_negative_strand, dBALF2_uninduced_positive_strand, dBALF2_uninduced_total_reads,regions=rta_regions)
rta_dBALF2_induced_coverage = extract_coverage(dBALF2_induced_negative_strand, dBALF2_induced_positive_strand, dBALF2_induced_total_reads,regions=rta_regions)
rta_dBALF2_induced_transcomp_coverage = extract_coverage(dBALF2_induced_transcomp_negative_strand, dBALF2_induced_transcomp_positive_strand, dBALF2_induced_transcomp_total_reads,regions=rta_regions)
rta_dBDLF4_induced_coverage = extract_coverage(dBDLF4_induced_negative_strand, dBDLF4_induced_positive_strand, dBDLF4_induced_total_reads,regions=rta_regions)
rta_dBDLF4_induced_transcomp_coverage = extract_coverage(dBDLF4_induced_transcomp_negative_strand, dBDLF4_induced_transcomp_positive_strand, dBDLF4_induced_transcomp_total_reads,regions=rta_regions)

plot_heatmap(rta_wt_induced_coverage, "rta_wt_induced_sense.denovo.CAGE.coverage", 'Reds', 8* (len(rta_regions)/len(denovo)))
plot_heatmap(rta_dOriLyt_induced_coverage, "rta_dOriLyt_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(rta_regions)/len(denovo)))
plot_heatmap(rta_dBALF2_uninduced_coverage, "rta_dBALF2_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(rta_regions)/len(denovo)))
plot_heatmap(rta_dBALF2_induced_coverage, "rta_dBALF2_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(rta_regions)/len(denovo)))
plot_heatmap(rta_dBALF2_induced_transcomp_coverage, "rta_dBALF2_induced_transcomp_sense.denovo.CAGE.coverage", 'Reds', 8*(len(rta_regions)/len(denovo)))
plot_heatmap(rta_dBDLF4_induced_coverage, "rta_dBDLF4_induced_sense.denovo.CAGE.coverage", "Reds", 8*(len(rta_regions)/len(denovo)))
plot_heatmap(rta_dBDLF4_induced_transcomp_coverage, "rta_dBDLF4_induced_transcomp_sense.denovo.CAGE.coverage", "Reds", 8*(len(rta_regions)/len(denovo)))


no_site_wt_induced_coverage = extract_coverage(wt_induced_negative_strand, wt_induced_positive_strand, wt_induced_total_reads,regions=no_site)
no_site_dOriLyt_induced_coverage = extract_coverage(dOriLyt_induced_negative_strand, dOriLyt_induced_positive_strand, dOriLyt_induced_total_reads,regions=no_site)
no_site_dBALF2_uninduced_coverage = extract_coverage(dBALF2_uninduced_negative_strand, dBALF2_uninduced_positive_strand, dBALF2_uninduced_total_reads,regions=no_site)
no_site_dBALF2_induced_coverage = extract_coverage(dBALF2_induced_negative_strand, dBALF2_induced_positive_strand, dBALF2_induced_total_reads,regions=no_site)
no_site_dBALF2_induced_transcomp_coverage = extract_coverage(dBALF2_induced_transcomp_negative_strand, dBALF2_induced_transcomp_positive_strand, dBALF2_induced_transcomp_total_reads,regions=no_site)
no_site_dBDLF4_induced_coverage = extract_coverage(dBDLF4_induced_negative_strand, dBDLF4_induced_positive_strand, dBDLF4_induced_total_reads,regions=no_site)
no_site_dBDLF4_induced_transcomp_coverage = extract_coverage(dBDLF4_induced_transcomp_negative_strand, dBDLF4_induced_transcomp_positive_strand, dBDLF4_induced_transcomp_total_reads,regions=no_site)

plot_heatmap(no_site_wt_induced_coverage, "no_site_wt_induced_sense.denovo.CAGE.coverage", 'Reds', 8* (len(no_site)/len(denovo)))
plot_heatmap(no_site_dOriLyt_induced_coverage, "no_site_dOriLyt_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(no_site)/len(denovo)))
plot_heatmap(no_site_dBALF2_uninduced_coverage, "no_site_dBALF2_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(no_site)/len(denovo)))
plot_heatmap(no_site_dBALF2_induced_coverage, "no_site_dBALF2_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(no_site)/len(denovo)))
plot_heatmap(no_site_dBALF2_induced_transcomp_coverage, "no_site_dBALF2_induced_transcomp_sense.denovo.CAGE.coverage", 'Reds', 8*(len(no_site)/len(denovo)))
plot_heatmap(no_site_dBDLF4_induced_coverage, "no_site_dBDLF4_induced_sense.denovo.CAGE.coverage", "Reds", 8*(len(no_site)/len(denovo)))
plot_heatmap(no_site_dBDLF4_induced_transcomp_coverage, "no_site_dBDLF4_induced_transcomp_sense.denovo.CAGE.coverage", "Reds", 8*(len(no_site)/len(denovo)))


