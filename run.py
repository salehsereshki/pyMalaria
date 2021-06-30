import methods as methods
import exon_boundry_plot as EBP
import gene_body_meth as GBM
import pie_plot as PP

NON_METH_TAG = 0.00000001

seq_address_pc = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PcynomolgiB_Genome.fasta'
meth_address_pc = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PC_merged.CX_report.txt'
annot_address_pc = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-50_PcynomolgiB.gff'
coverage_trsh = 3
organism_name = 'PC'

seq_address_pv = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PvivaxP01_Genome.fasta'
meth_address_pv = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PV_excluded23.CX_report.txt'
annot_address_pv = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PvivaxP01.gff'


seq_address_pf = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_Pfalciparum3D7_Genome.fasta'
meth_address_pf = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PF.CX_report.txt'
annot_address_pf = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_Pfalciparum3D7.gff'

seq_addresses = [seq_address_pc, seq_address_pv, seq_address_pf]
meth_addresses = [meth_address_pc, meth_address_pv, meth_address_pf]
annot_addresses = [annot_address_pc, annot_address_pv, annot_address_pf]
organism_names = ['PC', 'PV', 'PF']

for i in range(len(seq_addresses)):
    sequences = methods.readfasta(seq_addresses[i])
    methylations = methods.read_methylations(meth_addresses[i])
    # for PC gene methylation put coverage threshold 3 for PV put threshold 10
    if organism_names[i] == 'PC':
        coverage_thrshold = 3
    else:
        coverage_thrshold = 10
    meth_seq, context_seq = methods.make_meth_string(methylations, sequences, coverage_thrshold)
    annot_df = methods.read_annot(annot_addresses[i])
    genes_df = methods.subset_genes(annot_df)
    genes_seq = methods.make_gene_string(genes_df, sequences)
    EBP.plot_exon_boundry_3(organism_names[i], annot_df, meth_seq, boundry=150)
    EBP.plot_exon_boundry_5(organism_names[i], annot_df, meth_seq, boundry=150)
    GBM.plot_gene_body_meth(organism_names[i], meth_seq, genes_df, 5)
    PP.plot_meth_context_percentage(organism_names[i], methylations, coverage_thrshold)


