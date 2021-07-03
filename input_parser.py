from Bio import SeqIO
import pandas as pd
import numpy as np
import pyMalaria.constants as constants

def readfasta(address):
    cowpeaRecords = SeqIO.parse(address, "fasta")
    sequences = {}
    for chro in cowpeaRecords:
        sequences[chro.id] = chro.seq
    for i in sequences.keys():
        sequences[i] = sequences[i].upper()
    return sequences
def read_annot(address):
    annot_df = pd.read_table(address, sep='\t', comment='#')
    annot_df.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    return annot_df

def read_methylations(address):
    methylations = pd.read_table(address)
    methylations.columns = ['chr', 'position', 'strand', 'meth', 'unmeth', 'context', 'three']
    return methylations

def check_annot_chro_names(annot_df, sequences):
    annot_chrs = annot_df['chr'].unique()
    count = 0
    for chr in annot_chrs:
        if chr in sequences.keys():
            count += 1
    if len(annot_df['chr'].unique()) == count:
        return True
    else:
        return False

def subset_genes(annot_df):
    genes_df = annot_df[annot_df['type'] == 'gene']
    genes_df = genes_df.reset_index(drop=True)
    return genes_df[['chr', 'strand', 'start', 'end']]

def subset_exons(annot_df):
    gene_exon_df = annot_df[(annot_df['type'] == 'exon')]
    return gene_exon_df[['chr', 'strand', 'start', 'end']]

def make_gene_string(genes_df, sequences):
    genes_seq = {}
    for chr in sequences.keys():
        genes_seq[chr] = list(''.zfill(len(sequences[chr])))
    for index, row in genes_df.iterrows():
        gene_marker = '1'
        if row['strand'] == '-':
            gene_marker = '-1'
        for i in range(int(row['start'] - 1), int(row['end'] - 1)):
            genes_seq[row['chr']][i] = gene_marker
    return genes_seq

def make_exon_string(exons_df, sequences):
    exons_seq = {}
    for chr in sequences.keys():
        exons_seq[chr] = list(''.zfill(len(sequences[chr])))
    for index, row in exons_df.iterrows():
        gene_marker = '1'
        if row['strand'] == '-':
            gene_marker = '-1'
        for i in range(int(row['start']), int(row['end'])):
            exons_seq[row['chr']][i] = gene_marker
    return exons_seq


def make_meth_string(methylations, sequences, coverage_thrshld):
    methylations['mlevel'] = methylations['meth']/ (methylations['meth'] + methylations['unmeth'])
    methylations['coverage'] = methylations['meth'] + methylations['unmeth']
    methylations['mlevel'] = methylations['mlevel'].fillna(0)

    methylations.loc[(methylations.mlevel == 0),'mlevel'] = constants.NON_METH_TAG
    methylations.loc[(methylations.coverage < coverage_thrshld),'mlevel']= 0
    methylations.loc[(methylations.strand == '-'),'mlevel']= -1 * methylations.mlevel

    meth_seq = {}
    for chr in sequences.keys():
        meths = np.zeros(len(sequences[chr]))
        meth_subset = methylations[methylations['chr'] == chr]
        meths[[meth_subset['position'] - 1]] = meth_subset['mlevel']
        meth_seq[chr] = meths
    return meth_seq

def make_gene_plus_flanking_string(genes_df, sequences):
    make_gene_plus_flanking = {}
    for chr in sequences.keys():
        make_gene_plus_flanking[chr] = list(''.zfill(len(sequences[chr])))
    for index, row in genes_df.iterrows():
        gene_marker = '1'
        if row['strand'] == '-':
            gene_marker = '-1'
        for i in range(int(row['start']) - 500, int(row['end']) + 500):
            if i < len(make_gene_plus_flanking[row['chr']]) and i > 0:
                make_gene_plus_flanking[row['chr']][i] = gene_marker
    return make_gene_plus_flanking
