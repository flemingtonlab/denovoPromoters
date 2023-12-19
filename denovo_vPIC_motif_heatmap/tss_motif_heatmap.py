import sys
import pandas as pd
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import pybedtools
import numpy as np
import argparse


def generate_bedtool(df_path, upstream_bases, downstream_bases):

    df = pd.read_table(df_path, header=None)
    print(df.head())
    df['tss'] = (df[1] + df[2]) // 2
    left, right = [], []
    for strand, tss in zip(df[5], df['tss']):
        if strand == '+':
            left.append(tss - upstream_bases) 
            right.append(tss + downstream_bases)
        elif strand == '-':
            left.append(tss - downstream_bases)
            right.append(tss + upstream_bases)
    df['left'] = left
    df['right'] = right
    df = df[[0, 'left', 'right', 3, 4, 5]].sort_values(by=[0, 'left', 'right'])

    return df


def get_fa(bed_df, fa_path):
    '''Add the nucleotide sequence from each interval to the dataframe'''
    
    bed = pybedtools.BedTool.from_dataframe(bed_df) 
    bed = bed.sequence(fi=fa_path, s=True) # Fasta sequence id is '>chromosome:start-stop(strand)'
    bed_df.index = [f'{chrom}:{left}-{right}({strand})' for chrom, left, right, strand in zip(bed_df[0], bed_df['left'], bed_df['right'], bed_df[5])]
    seqs = [i.strip('>') for i in open(bed.seqfn).read().split('\n')]
    seqs_dict = {seqs[i]: seqs[i+1] for i in range(0, len(seqs)-1,2)}
    bed_df['sequence'] = bed_df.index.map(lambda x: seqs_dict[x])

    return bed_df


def get_motif_position(bed_df, motif): 

    bed_df['motifs'] = bed_df['sequence'].map(lambda x: [i.group() for i in motif.finditer(x)])
    bed_df['match_interval'] = bed_df['sequence'].map(lambda x: [i.span() for i in motif.finditer(x)])                  
    return bed_df


def build_heatmap_matrix(match_intervals, matrix, color=[255/255, 0, 0]):

    #expression = 1
    for ind, line in enumerate(match_intervals):
        for interval in line:
            matrix[ind, interval[0]] = color

    return matrix


def convert_match_interval(df, upstream_bases):
    mi = []
    for i in df['match_interval']:
        n = []
        for j in i:
            n.append((j[0] - upstream_bases, j[1] - upstream_bases))
        
        mi.append(n)
    df['match_interval'] = mi
    return df



denovo_bed_path = '/Users/nate/dnovo/beds/denovo_promoters/de_novo_promoters_Akata_BCR_plus_Mutu_Zta.bed'
canonical_bed_path = '/Users/nate/dnovo/TSS_from_mutu_ctl.bed'

upstream_bases = 100
downstream_bases = 100
motifs = [re.compile('TATT[TA]AA')]
genome_fasta_path = '/Users/nate/Documents/Genomes/Fastas/Genomic/hg38_plus_Akata_inverted.genome.fa'


range_length = upstream_bases + downstream_bases
for path in [denovo_bed_path, canonical_bed_path]:
    df = generate_bedtool(path, upstream_bases=upstream_bases, downstream_bases=downstream_bases)
    df = get_fa(df, fa_path=genome_fasta_path)
    length = upstream_bases + downstream_bases + 1
    matrix = np.ones([len(df.index), length, 3])
    background = [255/255, 245/255, 240/255]
    matrix [len(df.index) -1, length - 1] = background   

    for motif, color in zip(motifs,[[255/255, 0, 0]]):

        df = get_motif_position(df, motif)
        df = df.sort_values(4)
        matrix = build_heatmap_matrix(df['match_interval'], matrix, color)

    fig = plt.figure(figsize=(4,10))
    ax = plt.subplot()
    plt.imshow(matrix, aspect='auto', cmap='Reds') 
    ax.set_yticks([])
    ax.set_xlim(0, range_length-1)
    ax.set_xticks([0, upstream_bases, range_length])
    ax.set_xticklabels([f'-{upstream_bases}', 'TSS', f'+{downstream_bases}'])
    plt.savefig(f'{path}.svg', dpi=2000)
    plt.close('all')
    df = convert_match_interval(df, upstream_bases)
    df.to_csv(f'{path}.motif_heatmap.tsv', sep='\t')






