import methods
from os import listdir
from os.path import isfile, join
import numpy as np


pc_meth_adrss = [f for f in listdir('./') if isfile(join('./', f))]

## Coverage Thrshold is for calling a Cytosine methylated. A Cytosine only can be called methylated if the coverage is more than threshold.
## If a Cytosine has a coverage less than threshold and all the the reads at that loci are methylated it is still called unmethylated.
seq_address_pc = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PcynomolgiB_Genome.fasta'
meth_address_pc = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PC_merged.CX_report.txt'
annot_address_pc = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-50_PcynomolgiB.gff'
coverage_trsh = 3
organism_name = 'PC'

seq_address_pv = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PvivaxP01_Genome.fasta'
meth_address_pv = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PV_excluded23.CX_report.txt'
annot_address_pv = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_PvivaxP01.gff'
coverage_trsh = 10
organism_name = 'PV'

seq_address_pf = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_Pfalciparum3D7_Genome.fasta'
meth_address_pf = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PF.CX_report.txt'
annot_address_pf = '/Users/salehsereshki/PycharmProjects/pythonProject/malaria/PlasmoDB-48_Pfalciparum3D7.gff'
coverage_trsh = 10
organism_name = 'PF'

seq_addresses = [seq_address_pc, seq_address_pv, seq_address_pf]
meth_addresses = [meth_address_pc, meth_address_pv, meth_address_pf]
annot_addresses = [annot_address_pc, annot_address_pv, annot_address_pf]
organism_names = ['PC', 'PV', 'PF']

for i in range(len(seq_addresses)):
    sequences = methods.readfasta(seq_addresses[i])
    methylations = methods.read_methylations(meth_addresses[i])
    # for PC gene methylation put threshold 3 for PV put threshold 10
    if i == 0:
        meth_seq, context_seq = methods.make_meth_string(methylations, sequences, 3)
    else:
        meth_seq, context_seq = methods.make_meth_string(methylations, sequences, 10)
    annot_df = methods.read_annot(annot_addresses[i])
    genes_df = methods.subset_genes(annot_df)
    genes_seq = methods.make_gene_string(genes_df, sequences)

#get_gene_meth_count(seq_address2, meth_address2, annot_address2)

import methods as methods
sequences = methods.readfasta(seq_address_pv)
methylations = methods.read_methylations(meth_address_pv)
# for PC gene methylation put threshold 3 for PV put threshold 10
meth_seq, context_seq = methods.make_meth_string(methylations, sequences, coverage_trsh)
annot_df = methods.read_annot(annot_address_pv)
genes_df = methods.subset_genes(annot_df)
genes_seq = methods.make_gene_string(genes_df, sequences)

pc_chromosomes = ['DF157093', 'DF157094', 'DF157095', 'DF157096', 'DF157097', 'DF157098', 'DF157099', 'DF157100', 'DF157101', 'DF157102', 'DF157103', 'DF157104', 'DF157105', 'DF157106']
pv_chromosomes = ['PvP01_01_v1', 'PvP01_02_v1', 'PvP01_03_v1', 'PvP01_04_v1', 'PvP01_05_v1', 'PvP01_06_v1', 'PvP01_07_v1', 'PvP01_08_v1','PvP01_09_v1', 'PvP01_10_v1', 'PvP01_11_v1','PvP01_12_v1', 'PvP01_12_v1','PvP01_14_v1']
pf_chromosomes = ['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_04_v3', 'Pf3D7_05_v3', 'Pf3D7_06_v3', 'Pf3D7_07_v3', 'Pf3D7_08_v3', 'Pf3D7_09_v3', 'Pf3D7_10_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3', 'Pf3D7_13_v3', 'Pf3D7_14_v3']


import density_plot as DP
thresholds = [0.1]
for i in range(len(pv_chromosomes)):
        DP.plot_density_Cs(organism_name, pc_chromosomes[i], meth_seq, thresholds[0], i+1)
        for thrs in thresholds:
            DP.plot_density_methCs(organism_name, pv_chromosomes[i], meth_seq, thrs, i+1)


import exon_boundry_plot as EBP

EBP.plot_exon_boundry_5(organism_name, annot_df, meth_seq, boundry=150, from_file=False)
EBP.plot_exon_boundry_3(organism_name, annot_df, meth_seq, boundry=150, from_file=False)

import gene_body_meth as GBM
GBM.plot_gene_body_meth(organism_name, meth_seq, genes_df, 5)



import matplotlib.pyplot as plt
organisms = ['PF', 'PV', 'PC']
thresholds = [0.1]

for ogm in organisms:
    for thrs in thresholds:
        #res_meth_pct = np.load('/Users/salehsereshki/PycharmProjects/pythonProject/saved_data/'+ogm+'_thr_'+str(thrs)+'_comp_res_meth_pct.npy')
        res_pie_pct = np.load('/Users/salehsereshki/PycharmProjects/pythonProject/saved_data/'+ogm+'_thr_'+str(thrs)+'_comp_res_pie_pct.npy')
        plt.rcdefaults()
        font = {'size': 20}

        plt.rc('font', **font)
        labels = ['exons', 'introns', 'intergenic \n regions']
        explode = (0.05, 0.05, 0)
        fig1, ax1 = plt.subplots()
        colors=[rgb_to_hex((163, 65, 67)).upper(), rgb_to_hex((67, 94, 156)).upper(), rgb_to_hex((94, 150, 88)).upper()]
        colors = [colors[0], colors[2], colors[1]]
        res_pie_pct = [res_pie_pct[0], res_pie_pct[2], res_pie_pct[1]]
        labels = [labels[0], labels[2], labels[1]]
        ax1.pie(res_pie_pct, explode=explode, labels=None, autopct='%1.1f%%',
            shadow=False, rotatelabels = 270, textprops=dict(color="k"), colors=colors)
        #wedges, lbls, txts = ax1.pie(res_pie_pct, explode=explode, labels=labels, autopct='%1.1f%%',labeldistance = 1.1,
        #    shadow=False, startangle=0, rotatelabels = True, textprops=dict(color="k"), counterclock=True, colors=[rgb_to_hex((163, 65, 67)).upper(), rgb_to_hex((67, 94, 156)).upper(), rgb_to_hex((94, 150, 88)).upper()])

        #for label in lbls:
        #    label.set_horizontalalignment('left')

        ax1.legend(labels = labels,fontsize=16, ncol=3, loc = 'lower center', bbox_to_anchor=(0.5, -0.3) )

        ax1.axis('equal')
        plt.tight_layout()
        plt.savefig('compartment_meth_repartition' + str(ogm)+'_thr_'+str(thrs) + '.jpg', dpi=2000)
        plt.close()



