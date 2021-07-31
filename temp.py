import pyMalaria.configs as configs
import pyMalaria.gene_body_meth as GBM

config = configs.PF_config
organism_name, coverage_threshold, meth_threshold, chromosomes = configs.setup_config(config)
sequences, methylations, meth_seq, annot_df, genes_df = configs.initialize(config)

GBM.plot_gene_body_meth(organism_name, meth_seq, genes_df, 5)
