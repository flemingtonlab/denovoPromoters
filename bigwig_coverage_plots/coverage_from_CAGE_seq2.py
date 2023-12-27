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

akata_denovo_bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_all_Akata_BCR.bed'
mutu_denovo_bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_all_Mutu_Zta.bed'
canonical_bed_path = '/Users/nate/dnovo/beds/TSS_from_mutu_ctl2.bed'

uninduced_akata_paths_positive_strand = [f'/Users/nate/dnovo/cage/akata/TM{i}.hg38plusAkata_inverted-Signal.Unique.str1.out.wig.bw' for i in range(1,5)]
uninduced_akata_paths_negative_strand = [f'/Users/nate/dnovo/cage/akata/TM{i}.hg38plusAkata_inverted-Signal.Unique.str2.out.wig.minus.wig.bw' for i in range(1,5)]
induced_akata_paths_positive_strand = [f'/Users/nate/dnovo/cage/akata/TM{i}.hg38plusAkata_inverted-Signal.Unique.str1.out.wig.bw' for i in range(5,9)]
induced_akata_paths_negative_strand = [f'/Users/nate/dnovo/cage/akata/TM{i}.hg38plusAkata_inverted-Signal.Unique.str2.out.wig.minus.wig.bw' for i in range(5,9)]



ctl_mutu_paths_positive_strand = [f'/Users/nate/dnovo/cage/mutu/{i}.hg38plusAkata_inverted-Signal.Unique.str1.out.wig.bw' for i in ['MC1','MC2','MC4']]
ctl_mutu_paths_negative_strand = [f'/Users/nate/dnovo/cage/mutu/{i}.hg38plusAkata_inverted-Signal.Unique.str2.out.wig.minus.wig.bw' for i in ['MC1','MC2','MC4']]
zta_mutu_paths_positive_strand = [f'/Users/nate/dnovo/cage/mutu/{i}.hg38plusAkata_inverted-Signal.Unique.str1.out.wig.bw' for i in ['MZ1','MZ2','MZ4']]
zta_mutu_paths_negative_strand = [f'/Users/nate/dnovo/cage/mutu/{i}.hg38plusAkata_inverted-Signal.Unique.str2.out.wig.minus.wig.bw' for i in ['MZ1','MZ2','MZ4']]

ctl_mutu_mapped_reads = [16106692, 14336475, 15846688]
zta_mutu_mapped_reads = [11945092, 12183555, 12568417]

uninduced_akata_mapped_reads = [41083098,50524834,44814675,54712817]
induced_akata_mapped_reads= [12626134, 9045597,12798630, 29488729]

def import_bed(path, col_sort=4):
    regions = []
    with open(path) as bed_handle:
        for line in bed_handle:
            regions.append(line.strip('\n').split('\t'))
    random.shuffle(regions)
    regions.sort(key=lambda x: float(x[col_sort]))
    return regions


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
                vals = str2.values(region[0], start - upstream_bases, stop+downstream_bases)
            else:
                vals = str1.values(region[0], start - downstream_bases, stop+upstream_bases)
                vals = np.flip(vals)                        
            vals = np.abs(vals)
            vals = np.nan_to_num(np.array(vals), nan=0)
            m[ind] += (100000000 * vals) / total_coverage
        print(str1_path, "done")
    return m / len(strand1_paths)


def plot_heatmap(matrix, name, colormap='Reds',height=8):
    fig = plt.figure(figsize=(3,height))
    plt.imshow(matrix, cmap=colormap,aspect='auto', interpolation="gaussian", vmax=1)
    plt.xticks([])
    plt.yticks([])
    plt.savefig(name + '.heatmap.svg')
    plt.close('all') 

def plot_sumcurves(matrix1, name, color='r',cutoff_number=0, binsize=100, max_y=10):
    matrix1[matrix1 <= cutoff_number] = 0
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
    plt.xlim([0,len(y_vals)-1])
    plt.yticks([])
    plt.xticks([])
    print(np.max(sum_hist1[1:-1]) / np.shape(matrix1)[0],name)
    plt.savefig(name + '.summation.svg')

def get_zta_sites(denovo_regions):
    zta_vpic_df = pd.read_table('/Users/nate/Downloads/de_novo_promoter_paper/motifs_in_de_novo_promoters/4_vPIC_plus_Zta/4_TATTAAA_TATTTAA_and_Zta_ChIP_seq_Hammerschmidt/TATTAAA_plus_TATTTAA_and_Zta_ChIP_Hammerschmidt_output/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed_TATTAAA_plus_TATTTAA.bed.binary_motif_info.bed_Raji_Zta_ChIP_induced_pooledReps_summits.bed.no_EBV_plus_strand.bed.binary_motif_info.tsv', index_col=0)
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
        if i in l:
            two_plus_sites.append(i)
        l.append(i)

    no_site = [i for i in denovo_regions if i[3] not in all_three]
    vpic = [i for i in vpic if i[3] not in two_plus_sites]
    zta = [i for i in zta if i[3] not in two_plus_sites]
    rta = [i for i in rta if i[3] not in two_plus_sites]
    print(len(rta) + len(zta) + len(vpic) + len(no_site) + len(set(two_plus_sites)), len(denovo_regions))

    return rta, zta, vpic, no_site



mutu_denovo_regions = import_bed(mutu_denovo_bed_path)
akata_denovo_regions = import_bed(akata_denovo_bed_path)

mutu_rta, mutu_zta, mutu_vpic, mutu_no_site = get_sites(mutu_denovo_regions)
akata_rta, akata_zta, akata_vpic, akata_no_site = get_sites(akata_denovo_regions)



#####  AKATA  #####

vpic_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_vpic)
vpic_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_vpic)
vpic_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_vpic,antisense=True)
vpic_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_vpic,antisense=True)

zta_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_zta,antisense=False)
zta_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_zta,antisense=False)
zta_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_zta,antisense=True)
zta_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_zta,antisense=True)

rta_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_rta,antisense=True)
rta_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_rta,antisense=True)
rta_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_rta,antisense=False)
rta_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_rta,antisense=False)

no_site_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_no_site,antisense=True)
no_site_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_no_site,antisense=True)
no_site_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_no_site,antisense=False)
no_site_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_no_site,antisense=False)

plot_heatmap(vpic_akata_uninduced_coverage_sense, "vpic_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_vpic)/len(akata_denovo_regions)))
plot_heatmap(vpic_akata_induced_coverage_sense, "vpic_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_vpic)/len(akata_denovo_regions)))
plot_heatmap(vpic_akata_uninduced_coverage_antisense, "vpic_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_vpic)/len(akata_denovo_regions)))
plot_heatmap(vpic_akata_induced_coverage_antisense, "vpic_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_vpic)/len(akata_denovo_regions)))

plot_heatmap(zta_akata_uninduced_coverage_antisense, "zta_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_zta)/len(akata_denovo_regions)))
plot_heatmap(zta_akata_induced_coverage_antisense, "zta_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_zta)/len(akata_denovo_regions)))
plot_heatmap(zta_akata_uninduced_coverage_sense, "zta_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_zta)/len(akata_denovo_regions)))
plot_heatmap(zta_akata_induced_coverage_sense, "zta_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_zta)/len(akata_denovo_regions)))

plot_heatmap(rta_akata_uninduced_coverage_antisense, "rta_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_rta)/len(akata_denovo_regions)))
plot_heatmap(rta_akata_induced_coverage_antisense, "rta_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_rta)/len(akata_denovo_regions)))
plot_heatmap(rta_akata_uninduced_coverage_sense, "rta_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_rta)/len(akata_denovo_regions)))
plot_heatmap(rta_akata_induced_coverage_sense, "rta_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_rta)/len(akata_denovo_regions)))

plot_heatmap(no_site_akata_uninduced_coverage_antisense, "no_site_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_no_site)/len(akata_denovo_regions)))
plot_heatmap(no_site_akata_induced_coverage_antisense, "no_site_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_no_site)/len(akata_denovo_regions)))
plot_heatmap(no_site_akata_uninduced_coverage_sense, "no_site_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_no_site)/len(akata_denovo_regions)))
plot_heatmap(no_site_akata_induced_coverage_sense, "no_site_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_no_site)/len(akata_denovo_regions)))


plot_sumcurves(vpic_akata_uninduced_coverage_sense, "vpic_akata_uninduced_coverage_sense", 'r', max_y=5)
plot_sumcurves(vpic_akata_induced_coverage_sense, "vpic_akata_induced_coverage_sense", 'r', max_y=16)
plot_sumcurves(vpic_akata_uninduced_coverage_antisense, "vpic_akata_uninduced_coverage_antisense","b",max_y=5)
plot_sumcurves(vpic_akata_induced_coverage_antisense, "vpic_akata_induced_coverage_antisense", "b",max_y=16)

plot_sumcurves(zta_akata_uninduced_coverage_antisense, "zta_akata_uninduced_coverage_antisense", "b",max_y=5)
plot_sumcurves(zta_akata_induced_coverage_antisense, "zta_akata_induced_coverage_antisense", "b",max_y=16)
plot_sumcurves(zta_akata_uninduced_coverage_sense, "zta_akata_uninduced_coverage_sense", "r",max_y=5)
plot_sumcurves(zta_akata_induced_coverage_sense, "zta_akata_induced_coverage_sense", "r",max_y=16)

plot_sumcurves(rta_akata_uninduced_coverage_antisense,"rta_akata_uninduced_coverage_antisense", "b",max_y=5)
plot_sumcurves(rta_akata_induced_coverage_antisense, 'rta_akata_induced_coverage_antisense', 'b',max_y=16)
plot_sumcurves(rta_akata_uninduced_coverage_sense, 'rta_akata_uninduced_coverage_sense', 'r',max_y=5)
plot_sumcurves(rta_akata_induced_coverage_sense, 'rta_akata_induced_coverage_sense', 'r',max_y=16)

plot_sumcurves(no_site_akata_uninduced_coverage_antisense, "no_site_akata_uninduced_coverage_antisense", "b",max_y=5)
plot_sumcurves(no_site_akata_induced_coverage_antisense, "no_site_akata_induced_coverage_antisense", "b",max_y=16)
plot_sumcurves(no_site_akata_uninduced_coverage_sense, "no_site_akata_uninduced_coverage_sense", "r",max_y=5)
plot_sumcurves(no_site_akata_induced_coverage_sense, "no_site_akata_induced_coverage_sense", "r",max_y=16)


#####  MUTU  #####


vpic_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_vpic)
vpic_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand, zta_mutu_mapped_reads, regions=mutu_vpic)
vpic_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_vpic,antisense=True)
vpic_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_vpic,antisense=True)

zta_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_zta,antisense=True)
zta_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand, zta_mutu_mapped_reads, regions=mutu_zta,antisense=True)
zta_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_zta,antisense=False)
zta_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_zta,antisense=False)

rta_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_rta,antisense=True)
rta_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_rta,antisense=True)
rta_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_rta,antisense=False)
rta_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand, zta_mutu_mapped_reads, regions=mutu_rta,antisense=False)

no_site_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_no_site,antisense=True)
no_site_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_no_site,antisense=True)
no_site_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_no_site,antisense=False)
no_site_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand, zta_mutu_mapped_reads, regions=mutu_no_site,antisense=False)



plot_heatmap(vpic_mutu_ctl_coverage_sense, "vpic_mutu_ctl_coverage_sense", "Reds",8*(len(mutu_vpic)/len(mutu_denovo_regions)) )
plot_heatmap(vpic_mutu_zta_coverage_sense, "vpic_mutu_zta_coverage_sense", "Reds", 8*(len(mutu_vpic)/len(mutu_denovo_regions)))
plot_heatmap(vpic_mutu_ctl_coverage_antisense, "vpic_mutu_ctl_coverage_antisense","Blues",8*(len(mutu_vpic)/len(mutu_denovo_regions)))
plot_heatmap(vpic_mutu_zta_coverage_antisense, "vpic_mutu_zta_coverage_antisense", "Blues",8*(len(mutu_vpic)/len(mutu_denovo_regions)))

plot_heatmap(zta_mutu_ctl_coverage_antisense, "zta_mutu_ctl_coverage_antisense", "Blues",8*(len(mutu_zta)/len(mutu_denovo_regions)))
plot_heatmap(zta_mutu_zta_coverage_antisense, "zta_mutu_zta_coverage_antisense", "Blues",8*(len(mutu_zta)/len(mutu_denovo_regions)))
plot_heatmap(zta_mutu_ctl_coverage_sense, "zta_mutu_ctl_coverage_sense", "Reds",8*(len(mutu_zta)/len(mutu_denovo_regions)))
plot_heatmap(zta_mutu_zta_coverage_sense, "zta_mutu_zta_coverage_sense", "Reds",8*(len(mutu_zta)/len(mutu_denovo_regions)))

plot_heatmap(rta_mutu_ctl_coverage_antisense,"rta_mutu_ctl_coverage_antisense", "Blues",8*(len(mutu_rta)/len(mutu_denovo_regions)))
plot_heatmap(rta_mutu_zta_coverage_antisense, 'rta_mutu_zta_coverage_antisense', 'Blues',8*(len(mutu_rta)/len(mutu_denovo_regions)))
plot_heatmap(rta_mutu_ctl_coverage_sense, 'rta_mutu_ctl_coverage_sense', "Reds",8*(len(mutu_rta)/len(mutu_denovo_regions)))
plot_heatmap(rta_mutu_zta_coverage_sense, 'rta_mutu_zta_coverage_sense', "Reds",8*(len(mutu_rta)/len(mutu_denovo_regions)))

plot_heatmap(no_site_mutu_ctl_coverage_antisense, "no_site_mutu_ctl_coverage_antisense", "Blues",8*(len(mutu_no_site)/len(mutu_denovo_regions)))
plot_heatmap(no_site_mutu_zta_coverage_antisense, "no_site_mutu_zta_coverage_antisense", "Blues",8*(len(mutu_no_site)/len(mutu_denovo_regions)))
plot_heatmap(no_site_mutu_ctl_coverage_sense, "no_site_mutu_ctl_coverage_sense", "Reds",8*(len(mutu_no_site)/len(mutu_denovo_regions)))
plot_heatmap(no_site_mutu_zta_coverage_sense, "no_site_mutu_zta_coverage_sense", "Reds",8*(len(mutu_no_site)/len(mutu_denovo_regions)))


plot_sumcurves(vpic_mutu_ctl_coverage_sense, "vpic_mutu_ctl_coverage_sense", 'r', max_y=3,cutoff_number=0)
plot_sumcurves(vpic_mutu_zta_coverage_sense, "vpic_mutu_zta_coverage_sense", 'r', max_y=7,cutoff_number=0)
plot_sumcurves(vpic_mutu_ctl_coverage_antisense, "vpic_mutu_ctl_coverage_antisense","b",max_y=3,cutoff_number=0)
plot_sumcurves(vpic_mutu_zta_coverage_antisense, "vpic_mutu_zta_coverage_antisense", "b",max_y=7,cutoff_number=0)

plot_sumcurves(zta_mutu_ctl_coverage_antisense, "zta_mutu_ctl_coverage_antisense", "b",max_y=3,cutoff_number=0)
plot_sumcurves(zta_mutu_zta_coverage_antisense, "zta_mutu_zta_coverage_antisense", "b",max_y=7,cutoff_number=0)
plot_sumcurves(zta_mutu_ctl_coverage_sense, "zta_mutu_ctl_coverage_sense", "r",max_y=3,cutoff_number=0)
plot_sumcurves(zta_mutu_zta_coverage_sense, "zta_mutu_zta_coverage_sense", "r",max_y=7,cutoff_number=0)

plot_sumcurves(rta_mutu_ctl_coverage_antisense,"rta_mutu_ctl_coverage_antisense", "b",max_y=3,cutoff_number=0)
plot_sumcurves(rta_mutu_zta_coverage_antisense, 'rta_mutu_zta_coverage_antisense', 'b',max_y=7,cutoff_number=0)
plot_sumcurves(rta_mutu_ctl_coverage_sense, 'rta_mutu_ctl_coverage_sense', 'r',max_y=3,cutoff_number=0)
plot_sumcurves(rta_mutu_zta_coverage_sense, 'rta_mutu_zta_coverage_sense', 'r',max_y=7,cutoff_number=0)

plot_sumcurves(no_site_mutu_ctl_coverage_antisense, "no_site_mutu_ctl_coverage_antisense", "b",max_y=3,cutoff_number=0)
plot_sumcurves(no_site_mutu_zta_coverage_antisense, "no_site_mutu_zta_coverage_antisense", "b",max_y=7,cutoff_number=0)
plot_sumcurves(no_site_mutu_ctl_coverage_sense, "no_site_mutu_ctl_coverage_sense", "r",max_y=3,cutoff_number=0)
plot_sumcurves(no_site_mutu_zta_coverage_sense, "no_site_mutu_zta_coverage_sense", "r",max_y=7,cutoff_number=0)



####################    SHIFTED

mutu_denovo_shifted =[]
akata_denovo_shifted = []
for i in mutu_denovo_regions:
    if i[5] == '-':
        mutu_denovo_shifted.append([i[0], int(i[1])+5000, int(i[2])+5000, i[3], i[4], i[5]])
    else:
        mutu_denovo_shifted.append([i[0], int(i[1])-5000, int(i[2])-5000, i[3], i[4], i[5]])

for i in akata_denovo_regions:
    if i[5] == '-':
        akata_denovo_shifted.append([i[0], int(i[1])+5000, int(i[2])+5000, i[3], i[4], i[5]])
    else: 
        akata_denovo_shifted.append([i[0], int(i[1])-5000, int(i[2])-5000, i[3], i[4], i[5]])


mutu_rta_shifted, mutu_zta_shifted, mutu_vpic_shifted, mutu_no_site_shifted = get_sites(mutu_denovo_shifted)
akata_rta_shifted, akata_zta_shifted, akata_vpic_shifted, akata_no_site_shifted = get_sites(akata_denovo_shifted)


shifted_vpic_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_vpic_shifted)
shifted_vpic_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_vpic_shifted)
shifted_vpic_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand,ctl_mutu_mapped_reads, regions=mutu_vpic_shifted,antisense=True)
shifted_vpic_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_vpic_shifted,antisense=True)

shifted_zta_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_zta_shifted,antisense=True)
shifted_zta_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_zta_shifted,antisense=True)
shifted_zta_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_zta_shifted,antisense=False)
shifted_zta_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand, zta_mutu_mapped_reads, regions=mutu_zta_shifted,antisense=False)

shifted_rta_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_rta_shifted,antisense=True)
shifted_rta_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_rta_shifted,antisense=True)
shifted_rta_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_rta_shifted,antisense=False)
shifted_rta_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_rta_shifted,antisense=False)

shifted_no_site_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_no_site_shifted,antisense=True)
shifted_no_site_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_no_site_shifted,antisense=True)
shifted_no_site_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=mutu_no_site_shifted,antisense=False)
shifted_no_site_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand,  zta_mutu_mapped_reads, regions=mutu_no_site_shifted,antisense=False)


plot_heatmap(shifted_vpic_mutu_ctl_coverage_sense, "shifted_vpic_mutu_ctl_coverage_sense", "Reds",8*(len(mutu_vpic_shifted)/len(mutu_denovo_shifted)) )
plot_heatmap(shifted_vpic_mutu_zta_coverage_sense, "shifted_vpic_mutu_zta_coverage_sense", "Reds", 8*(len(mutu_vpic_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_vpic_mutu_ctl_coverage_antisense, "shifted_vpic_mutu_ctl_coverage_antisense","Blues",8*(len(mutu_vpic_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_vpic_mutu_zta_coverage_antisense, "shifted_vpic_mutu_zta_coverage_antisense", "Blues",8*(len(mutu_vpic_shifted)/len(mutu_denovo_shifted)))

plot_heatmap(shifted_zta_mutu_ctl_coverage_antisense, "shifted_zta_mutu_ctl_coverage_antisense", "Blues",8*(len(mutu_zta_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_zta_mutu_zta_coverage_antisense, "shifted_zta_mutu_zta_coverage_antisense", "Blues",8*(len(mutu_zta_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_zta_mutu_ctl_coverage_sense, "shifted_zta_mutu_ctl_coverage_sense", "Reds",8*(len(mutu_zta_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_zta_mutu_zta_coverage_sense, "shifted_zta_mutu_zta_coverage_sense", "Reds",8*(len(mutu_zta_shifted)/len(mutu_denovo_shifted)))

plot_heatmap(shifted_rta_mutu_ctl_coverage_antisense,"shifted_rta_mutu_ctl_coverage_antisense", "Blues",8*(len(mutu_rta_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_rta_mutu_zta_coverage_antisense, 'shifted_rta_mutu_zta_coverage_antisense', 'Blues',8*(len(mutu_rta_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_rta_mutu_ctl_coverage_sense, 'shifted_rta_mutu_ctl_coverage_sense', "Reds",8*(len(mutu_rta_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_rta_mutu_zta_coverage_sense, 'shifted_rta_mutu_zta_coverage_sense', "Reds",8*(len(mutu_rta_shifted)/len(mutu_denovo_shifted)))

plot_heatmap(shifted_no_site_mutu_ctl_coverage_antisense, "shifted_no_site_mutu_ctl_coverage_antisense", "Blues",8*(len(mutu_no_site_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_no_site_mutu_zta_coverage_antisense, "shifted_no_site_mutu_zta_coverage_antisense", "Blues",8*(len(mutu_no_site_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_no_site_mutu_ctl_coverage_sense, "shifted_no_site_mutu_ctl_coverage_sense", "Reds",8*(len(mutu_no_site_shifted)/len(mutu_denovo_shifted)))
plot_heatmap(shifted_no_site_mutu_zta_coverage_sense, "shifted_no_site_mutu_zta_coverage_sense", "Reds",8*(len(mutu_no_site_shifted)/len(mutu_denovo_shifted)))


plot_sumcurves(shifted_vpic_mutu_ctl_coverage_sense, "shifted_vpic_mutu_ctl_coverage_sense", 'r', max_y=3,cutoff_number=0)
plot_sumcurves(shifted_vpic_mutu_zta_coverage_sense, "shifted_vpic_mutu_zta_coverage_sense", 'r', max_y=7,cutoff_number=0)
plot_sumcurves(shifted_vpic_mutu_ctl_coverage_antisense, "shifted_vpic_mutu_ctl_coverage_antisense","b",max_y=3,cutoff_number=0)
plot_sumcurves(shifted_vpic_mutu_zta_coverage_antisense, "shifted_vpic_mutu_zta_coverage_antisense", "b",max_y=7,cutoff_number=0)

plot_sumcurves(shifted_zta_mutu_ctl_coverage_antisense, "shifted_zta_mutu_ctl_coverage_antisense", "b",max_y=3,cutoff_number=0)
plot_sumcurves(shifted_zta_mutu_zta_coverage_antisense, "shifted_zta_mutu_zta_coverage_antisense", "b",max_y=7,cutoff_number=0)
plot_sumcurves(shifted_zta_mutu_ctl_coverage_sense, "shifted_zta_mutu_ctl_coverage_sense", "r",max_y=3,cutoff_number=0)
plot_sumcurves(shifted_zta_mutu_zta_coverage_sense, "shifted_zta_mutu_zta_coverage_sense", "r",max_y=7,cutoff_number=0)

plot_sumcurves(shifted_rta_mutu_ctl_coverage_antisense,"shifted_rta_mutu_ctl_coverage_antisense", "b",max_y=3,cutoff_number=0)
plot_sumcurves(shifted_rta_mutu_zta_coverage_antisense, 'shifted_rta_mutu_zta_coverage_antisense', 'b',max_y=7,cutoff_number=0)
plot_sumcurves(shifted_rta_mutu_ctl_coverage_sense, 'shifted_rta_mutu_ctl_coverage_sense', 'r',max_y=3,cutoff_number=0)
plot_sumcurves(shifted_rta_mutu_zta_coverage_sense, 'shifted_rta_mutu_zta_coverage_sense', 'r',max_y=7,cutoff_number=0)

plot_sumcurves(shifted_no_site_mutu_ctl_coverage_antisense, "shifted_no_site_mutu_ctl_coverage_antisense", "b",max_y=3,cutoff_number=0)
plot_sumcurves(shifted_no_site_mutu_zta_coverage_antisense, "shifted_no_site_mutu_zta_coverage_antisense", "b",max_y=7,cutoff_number=0)
plot_sumcurves(shifted_no_site_mutu_ctl_coverage_sense, "shifted_no_site_mutu_ctl_coverage_sense", "r",max_y=3,cutoff_number=0)
plot_sumcurves(shifted_no_site_mutu_zta_coverage_sense, "shifted_no_site_mutu_zta_coverage_sense", "r",max_y=7,cutoff_number=0)



########### SHIFTED AKATA

shifted_vpic_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_vpic_shifted)
shifted_vpic_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_vpic_shifted)
shifted_vpic_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_vpic_shifted,antisense=True)
shifted_vpic_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand, induced_akata_mapped_reads, regions=akata_vpic_shifted,antisense=True)

shifted_zta_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_zta_shifted,antisense=True)
shifted_zta_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand, induced_akata_mapped_reads, regions=akata_zta_shifted,antisense=True)
shifted_zta_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_zta_shifted,antisense=False)
shifted_zta_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_zta_shifted,antisense=False)

shifted_rta_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_rta_shifted,antisense=True)
shifted_rta_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand, induced_akata_mapped_reads, regions=akata_rta_shifted,antisense=True)
shifted_rta_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_rta_shifted,antisense=False)
shifted_rta_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,induced_akata_mapped_reads, regions=akata_rta_shifted,antisense=False)

shifted_no_site_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_no_site_shifted,antisense=True)
shifted_no_site_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand,  induced_akata_mapped_reads, regions=akata_no_site_shifted,antisense=True)
shifted_no_site_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=akata_no_site_shifted,antisense=False)
shifted_no_site_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand, induced_akata_mapped_reads, regions=akata_no_site_shifted,antisense=False)


plot_heatmap(shifted_vpic_akata_uninduced_coverage_sense, "shifted_vpic_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8* (len(akata_vpic_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_vpic_akata_induced_coverage_sense, "shifted_vpic_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_vpic_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_vpic_akata_uninduced_coverage_antisense, "shifted_vpic_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_vpic_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_vpic_akata_induced_coverage_antisense, "shifted_vpic_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_vpic_shifted)/len(akata_denovo_shifted)))

plot_heatmap(shifted_zta_akata_uninduced_coverage_antisense, "shifted_zta_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_zta_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_zta_akata_induced_coverage_antisense, "shifted_zta_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_zta_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_zta_akata_uninduced_coverage_sense, "shifted_zta_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_zta_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_zta_akata_induced_coverage_sense, "shifted_zta_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_zta_shifted)/len(akata_denovo_shifted)))

plot_heatmap(shifted_rta_akata_uninduced_coverage_antisense, "shifted_rta_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_rta_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_rta_akata_induced_coverage_antisense, "shifted_rta_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_rta_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_rta_akata_uninduced_coverage_sense, "shifted_rta_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_rta_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_rta_akata_induced_coverage_sense, "shifted_rta_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_rta_shifted)/len(akata_denovo_shifted)))

plot_heatmap(shifted_no_site_akata_uninduced_coverage_antisense, "shifted_no_site_akata_uninduced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_no_site_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_no_site_akata_induced_coverage_antisense, "shifted_no_site_akata_induced_antisense.denovo.CAGE.coverage", 'Blues', 8*(len(akata_no_site_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_no_site_akata_uninduced_coverage_sense, "shifted_no_site_akata_uninduced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_no_site_shifted)/len(akata_denovo_shifted)))
plot_heatmap(shifted_no_site_akata_induced_coverage_sense, "shifted_no_site_akata_induced_sense.denovo.CAGE.coverage", 'Reds', 8*(len(akata_no_site_shifted)/len(akata_denovo_shifted)))

plot_sumcurves(shifted_vpic_akata_uninduced_coverage_sense, "shifted_vpic_akata_uninduced_coverage_sense", 'r', max_y=5)
plot_sumcurves(shifted_vpic_akata_induced_coverage_sense, "shifted_vpic_akata_induced_coverage_sense", 'r', max_y=16)
plot_sumcurves(shifted_vpic_akata_uninduced_coverage_antisense, "shifted_vpic_akata_uninduced_coverage_antisense","b",max_y=5)
plot_sumcurves(shifted_vpic_akata_induced_coverage_antisense, "shifted_vpic_akata_induced_coverage_antisense", "b",max_y=16)

plot_sumcurves(shifted_zta_akata_uninduced_coverage_antisense, "shifted_zta_akata_uninduced_coverage_antisense", "b",max_y=5)
plot_sumcurves(shifted_zta_akata_induced_coverage_antisense, "shifted_zta_akata_induced_coverage_antisense", "b",max_y=16)
plot_sumcurves(shifted_zta_akata_uninduced_coverage_sense, "shifted_zta_akata_uninduced_coverage_sense", "r",max_y=5)
plot_sumcurves(shifted_zta_akata_induced_coverage_sense, "shifted_zta_akata_induced_coverage_sense", "r",max_y=16)

plot_sumcurves(shifted_rta_akata_uninduced_coverage_antisense,"shifted_rta_akata_uninduced_coverage_antisense", "b",max_y=5)
plot_sumcurves(shifted_rta_akata_induced_coverage_antisense, 'shifted_rta_akata_induced_coverage_antisense', 'b',max_y=16)
plot_sumcurves(shifted_rta_akata_uninduced_coverage_sense, 'shifted_rta_akata_uninduced_coverage_sense', 'r',max_y=5)
plot_sumcurves(shifted_rta_akata_induced_coverage_sense, 'shifted_rta_akata_induced_coverage_sense', 'r',max_y=16)

plot_sumcurves(shifted_no_site_akata_uninduced_coverage_antisense, "shifted_no_site_akata_uninduced_coverage_antisense", "b",max_y=5)
plot_sumcurves(shifted_no_site_akata_induced_coverage_antisense, "shifted_no_site_akata_induced_coverage_antisense", "b",max_y=16)
plot_sumcurves(shifted_no_site_akata_uninduced_coverage_sense, "shifted_no_site_akata_uninduced_coverage_sense", "r",max_y=5)
plot_sumcurves(shifted_no_site_akata_induced_coverage_sense, "shifted_no_site_akata_induced_coverage_sense", "r",max_y=16)



canonical_bed_path = '/Users/nate/dnovo/beds/TSS_from_mutu_ctl2.bed'
canonical = import_bed(canonical_bed_path,col_sort=6)

canonical_mutu_ctl_coverage_sense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=canonical)
canonical_mutu_ctl_coverage_antisense = extract_coverage(ctl_mutu_paths_negative_strand, ctl_mutu_paths_positive_strand, ctl_mutu_mapped_reads, regions=canonical,antisense=True)

canonical_mutu_zta_coverage_sense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand, zta_mutu_mapped_reads, regions=canonical)
canonical_mutu_zta_coverage_antisense = extract_coverage(zta_mutu_paths_negative_strand, zta_mutu_paths_positive_strand, zta_mutu_mapped_reads, regions=canonical,antisense=True)

canonical_akata_uninduced_coverage_sense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=canonical)
canonical_akata_uninduced_coverage_antisense = extract_coverage(uninduced_akata_paths_negative_strand, uninduced_akata_paths_positive_strand, uninduced_akata_mapped_reads, regions=canonical,antisense=True)

canonical_akata_induced_coverage_sense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand, induced_akata_mapped_reads, regions=canonical)
canonical_akata_induced_coverage_antisense = extract_coverage(induced_akata_paths_negative_strand, induced_akata_paths_positive_strand, induced_akata_mapped_reads, regions=canonical,antisense=True)


plot_heatmap(canonical_mutu_ctl_coverage_sense, "canonical_mutu_ctl_coverage_sense.denovo.CAGE.coverage", 'Reds', 10)
plot_heatmap(canonical_mutu_ctl_coverage_antisense, "canonical_mutu_ctl_coverage_antisense.denovo.CAGE.coverage", 'Blues', 10)
plot_heatmap(canonical_mutu_zta_coverage_sense, "canonical_mutu_zta_coverage_sense.denovo.CAGE.coverage", 'Reds', 10)
plot_heatmap(canonical_mutu_zta_coverage_antisense, "canonical_mutu_zta_coverage_antisense.denovo.CAGE.coverage", 'Blues', 10)


plot_heatmap(canonical_akata_uninduced_coverage_sense, "canonical_akata_uninduced_coverage_sense.denovo.CAGE.coverage", 'Reds', 10)
plot_heatmap(canonical_akata_uninduced_coverage_antisense, "canonical_akata_uninduced_coverage_antisense.denovo.CAGE.coverage", 'Blues', 10)
plot_heatmap(canonical_akata_induced_coverage_sense, "canonical_akata_induced_coverage_sense.denovo.CAGE.coverage", 'Reds', 10)
plot_heatmap(canonical_akata_induced_coverage_antisense, "canonical_akata_induced_coverage_antisense.denovo.CAGE.coverage", 'Blues', 10)



plot_sumcurves(canonical_mutu_ctl_coverage_sense, "canonical_mutu_ctl_coverage_sense.denovo.CAGE.coverage", 'r', max_y=4)
plot_sumcurves(canonical_mutu_ctl_coverage_antisense, "canonical_mutu_ctl_coverage_antisense.denovo.CAGE.coverage", 'b', max_y=4)
plot_sumcurves(canonical_mutu_zta_coverage_sense, "canonical_mutu_zta_coverage_sense.denovo.CAGE.coverage", 'r', max_y=4)
plot_sumcurves(canonical_mutu_zta_coverage_antisense, "canonical_mutu_zta_coverage_antisense.denovo.CAGE.coverage", 'b', max_y=4)


plot_sumcurves(canonical_akata_uninduced_coverage_sense, "canonical_akata_uninduced_coverage_sense.denovo.CAGE.coverage", 'r', max_y=4)
plot_sumcurves(canonical_akata_uninduced_coverage_antisense, "canonical_akata_uninduced_coverage_antisense.denovo.CAGE.coverage", 'b', max_y=4)
plot_sumcurves(canonical_akata_induced_coverage_sense, "canonical_akata_induced_coverage_sense.denovo.CAGE.coverage", 'r', max_y=4)
plot_sumcurves(canonical_akata_induced_coverage_antisense, "canonical_akata_induced_coverage_antisense.denovo.CAGE.coverage", 'b', max_y=4)

