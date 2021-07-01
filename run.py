import methods as methods
import exon_boundry_plot as EBP
import gene_body_meth as GBM
import pie_plot as PP

import pyMalaria.temp as temp

config = temp.PV_config

config['meth_threshold'] = 0.1

organism_name = config['organism_name']
coverage_threshold = config['coverage_threshold']
meth_threshold = config['meth_threshold']
chromosomes = config['chromosomes']


sequences, methylations, meth_seq, annot_df, genes_df = temp.initialize(config)

import pyMalaria.find_high_meth_regions as FHMR



FHMR.run(meth_seq, annot_df, organism_name, chromosomes, meth_threshold = 0.1, meth_C_in_kbp = 30)

