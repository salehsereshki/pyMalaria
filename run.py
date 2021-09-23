import methods as methods
import exon_boundry_plot as EBP
import gene_body_meth as GBM
import pie_plot as PP
import constants as constants
import configs as configs
import input_parser as IP
import matplotlib.pyplot as plt
import numpy as np

config = configs.PV_config
organism_name, coverage_threshold, meth_threshold, chromosomes = configs.setup_config(config)
sequences, methylations, meth_seq, annot_df, genes_df = configs.initialize(config)

meth_count_seq, unmeth_count_seq = IP.make_meth_count_string(methylations, sequences)

import find_high_meth_regions as FHMR



FHMR.run(meth_seq, annot_df, organism_name, chromosomes, meth_threshold = 0.1, meth_C_in_kbp = 30)


import gene_body_meth as GBM

import temp as tmp

tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 1, threshold = 0.1)
tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 2, threshold = 0.1)
tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 3, threshold = 0.1)
tmp.plot_gene_body_meth(organism_name, meth_seq, genes_df, 5, 4, threshold=0.1)

tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 5, threshold = 0.1)


FSs = [500]
thrshlds = [0.01, 0.05, 0.1, 0.15, 0.2]
config_list = [configs.PF_config, configs.PV_config, configs.PC_config]

for config in config_list:
    organism_name, coverage_threshold, meth_threshold, chromosomes = configs.setup_config(config)
    sequences, methylations, meth_seq, annot_df, genes_df = configs.initialize(config)
    for fs in FSs:
        for thrs in thrshlds:
            tmp.plot_gene_body_meth(organism_name, meth_seq, genes_df, 5, 4, threshold=thrs, flanking_size=fs)


FSs = [250, 500, 1000, 1500, 2000]
config_list = [configs.PF_config, configs.PV_config, configs.PC_config]

for config in config_list:
    organism_name, coverage_threshold, meth_threshold, chromosomes = configs.setup_config(config)
    sequences, methylations, meth_seq, annot_df, genes_df = configs.initialize(config)
    meth_count_seq, unmeth_count_seq = IP.make_meth_count_string(methylations, sequences)
    for fs in FSs:
        tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 5, threshold=0.1, flanking_size=fs)

import pyMalaria.expression_analysis as EA
exp_filters = [0.1, 0.25, 0.33]
thrs = 0.1
fs = 500
for exp_prcnt in exp_filters:
    exp_df = EA.get_most_expressed_genes(annot_df, configs.PV_expression, exp_prcnt)
    tmp.plot_gene_body_meth(organism_name, meth_seq, exp_df, 5, 4, threshold=thrs, flanking_size=fs, exp_filter=exp_prcnt)
    exp_df = EA.get_most_expressed_genes(annot_df, configs.PV_expression, -1 * exp_prcnt)
    tmp.plot_gene_body_meth(organism_name, meth_seq, exp_df, 5, 4, threshold=thrs, flanking_size=fs, exp_filter=-1*exp_prcnt)

import pyMalaria.gbm_computations as gbmc
import pyMalaria.expression_analysis as EA
exp_filters = [0.1, 0.25, 0.33]
thrs = 0.1
fs = 500
for exp_prcnt in exp_filters:
    exp_df = EA.get_most_expressed_genes(annot_df, configs.PV_expression, exp_prcnt)
    gbmc.get_gene_meth(meth_seq, exp_df,  5, exp_prcnt, threshold=0.1, flanking_size=500)
    exp_df = EA.get_most_expressed_genes(annot_df, configs.PV_expression, -1 * exp_prcnt)
    gbmc.get_gene_meth(meth_seq, exp_df,  5, -1 * exp_prcnt, threshold=0.1, flanking_size=500)


def checklen(df):
    for i in range(len(df)):
        if len(df[i]) != 15:
            print(i)


def gbm_df_plt(df, file_name):
    final_p = np.nanmean(df, axis=0)
    plt.ylim([0, 0.015])
    plt.plot(range(0,len(final_p)), final_p, color='blue', linewidth=4.0)
    my_path = '/Users/salehsereshki/Desktop/Malaria/me_high_low_dfs/coverage5/dfs/'
    plt.savefig(my_path + 'genebody_' + file_name + '.jpg', dpi=2000)
    plt.close()
from os import listdir
from os.path import isfile, join
my_path = '/Users/salehsereshki/Desktop/Malaria/me_high_low_dfs/coverage5/dfs/'
onlyfiles = [f for f in listdir(my_path) if isfile(join(my_path, f))]
for name in onlyfiles:
    if 'df' in name:
        df = np.load(my_path + name,  allow_pickle=True)
        gbm_df_plt(df, name)
