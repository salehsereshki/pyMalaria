import pyMalaria.input_parser as input_parser

import pyMalaria.density_plot as DP
import pyMalaria.genome_area_dist_plot as GADP
import pyMalaria.exon_boundry_plot as EBP
import pyMalaria.gc_content as GCC
import pyMalaria.gene_body_meth as GBM
import pyMalaria.pie_plot as PP

pc_chromosomes = ['DF157093', 'DF157094', 'DF157095', 'DF157096', 'DF157097', 'DF157098', 'DF157099', 'DF157100', 'DF157101', 'DF157102', 'DF157103', 'DF157104', 'DF157105', 'DF157106']
pv_chromosomes = ['PvP01_01_v1', 'PvP01_02_v1', 'PvP01_03_v1', 'PvP01_04_v1', 'PvP01_05_v1', 'PvP01_06_v1', 'PvP01_07_v1', 'PvP01_08_v1','PvP01_09_v1', 'PvP01_10_v1', 'PvP01_11_v1','PvP01_12_v1', 'PvP01_12_v1','PvP01_14_v1']
pf_chromosomes = ['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_04_v3', 'Pf3D7_05_v3', 'Pf3D7_06_v3', 'Pf3D7_07_v3', 'Pf3D7_08_v3', 'Pf3D7_09_v3', 'Pf3D7_10_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3', 'Pf3D7_13_v3', 'Pf3D7_14_v3']


PC_config = {
    'seq_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PcynomolgiB_Genome.fasta',
    'meth_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PC_merged.CX_report.txt',
    'annot_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-50_PcynomolgiB.gff',
    'coverage_threshold': 3,
    'organism_name': 'PC',
    'chromosomes': pc_chromosomes
}

PV_config = {
    'seq_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PvivaxP01_Genome.fasta',
    'meth_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PV_excluded23.CX_report.txt',
    'annot_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PvivaxP01.gff',
    'coverage_threshold': 10,
    'organism_name': 'PV',
    'chromosomes': pv_chromosomes
}

PF_config = {
    'seq_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_Pfalciparum3D7_Genome.fasta',
    'meth_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PF.CX_report.txt',
    'annot_address': '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_Pfalciparum3D7.gff',
    'coverage_threshold': 10,
    'organism_name': 'PF',
    'chromosomes': pf_chromosomes
}

PV_expression = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PvivaxP01_rnaseq.txt'


def setup_config(config):
    config['meth_threshold'] = 0.1

    organism_name = config['organism_name']
    coverage_threshold = config['coverage_threshold']
    meth_threshold = config['meth_threshold']
    chromosomes = config['chromosomes']
    return organism_name, coverage_threshold, meth_threshold, chromosomes



def initialize(config):
    sequences = input_parser.readfasta(config['seq_address'])
    methylations = input_parser.read_methylations(config['meth_address'])
    meth_seq = input_parser.make_meth_string(methylations, sequences, config['coverage_threshold'])
    annot_df = input_parser.read_annot(config['annot_address'])
    genes_df = input_parser.subset_genes(annot_df)
    return sequences, methylations, meth_seq, annot_df, genes_df




def run(config):
    organism_name, coverage_threshold, meth_threshold, chromosomes = setup_config(config)
    sequences, methylations, meth_seq, annot_df, genes_df = initialize(config)
    # DP.make_all_chromosome_density_plot(chromosomes, organism_name, meth_seq)

    GADP.plot_meth_percent_genome_areas(organism_name, annot_df, sequences, methylations, coverage_threshold,
                                thrshold=meth_threshold, from_file=False)

    EBP.plot_exon_boundry_5(organism_name, annot_df, meth_seq, boundry=150, from_file=False)
    EBP.plot_exon_boundry_3(organism_name, annot_df, meth_seq, boundry=150, from_file=False)

    GCC.plot_cg_content_percentage(sequences, genes_df, 5, organism_name)

    GBM.plot_gene_body_meth(organism_name, meth_seq, genes_df, 5)

    PP.plot_meth_context_percentage(organism_name, methylations, coverage_threshold)
