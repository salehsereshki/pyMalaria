import methods as methods
import exon_boundry_plot as EBP
import gene_body_meth as GBM
import pie_plot as PP
import constants as constants
import pyMalaria.configs as configs
import pyMalaria.input_parser as IP

config = configs.PC_config
organism_name, coverage_threshold, meth_threshold, chromosomes = configs.setup_config(config)
sequences, methylations, meth_seq, annot_df, genes_df = configs.initialize(config)

meth_count_seq, unmeth_count_seq = IP.make_meth_count_string(methylations, sequences)

import pyMalaria.find_high_meth_regions as FHMR



FHMR.run(meth_seq, annot_df, organism_name, chromosomes, meth_threshold = 0.1, meth_C_in_kbp = 30)


import pyMalaria.gene_body_meth as GBM

import pyMalaria.temp as tmp

tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 1, threshold = 0.1)
tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 2, threshold = 0.1)
tmp.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, 3, threshold = 0.1)
tmp.plot_gene_body_meth(organism_name, meth_seq, genes_df, 5, 4, threshold = 0.1)

GBM.plot_gene_body_meth(organism_name, [meth_count_seq, unmeth_count_seq], genes_df, 5, threshold = 0.1)
