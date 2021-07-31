import methods as methods
import exon_boundry_plot as EBP
import gene_body_meth as GBM
import pie_plot as PP
import constants as constants
import pyMalaria.configs as configs

config = configs.PV_config
organism_name, coverage_threshold, meth_threshold, chromosomes = configs.setup_config(config)
sequences, methylations, meth_seq, annot_df, genes_df = configs.initialize(config)


import pyMalaria.find_high_meth_regions as FHMR



FHMR.run(meth_seq, annot_df, organism_name, chromosomes, meth_threshold = 0.1, meth_C_in_kbp = 30)


import pyMalaria.gene_body_meth as GBM

GBM.plot_gene_body_meth(organism_name, meth_seq, genes_df, 5, threshold = 0.1)
