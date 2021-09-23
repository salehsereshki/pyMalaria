import configs
import methods
import numpy as np
from scipy.stats import pearsonr
import input_parser as input_parser
import constants as constants
import glob

def get_pearson_mtx(meht_address_lst, seq_address, config):
    sample_meths = []
    sequences = input_parser.readfasta(seq_address)
    chros = list(sequences.keys())
    for i in meht_address_lst:
        methylations = input_parser.read_methylations(meht_address_lst[i])
        meth_seq = input_parser.make_meth_string_plus_no_coverage(methylations, sequences, config['coverage_threshold'])
        meth_str = np.asarray([])
        for chro in (chros):
            meth_str = np.concatenate((meth_str, meth_seq[chro]), axis=0)
        sample_meths.append(meth_str)

        print('smple meth seq generated', i)
    print('sample meth seq generating finished')
    for i in range(len(sample_meths)):
        sample_meths[i] = np.asarray(sample_meths[i], dtype ='float64')
        sample_meths[i] = sample_meths[i][sample_meths[i] != 0]
        sample_meths[i] = np.abs(sample_meths[i])
        sample_meths[i] = np.where(sample_meths[i] == constants.NON_METH_TAG, 0, sample_meths[i])
    print('meth sample correction finished')
    corr_mtx = []
    for i in range(len(sample_meths)):
        row = []
        for j in range(len(sample_meths)):
            corr, _ = pearsonr(sample_meths[i], sample_meths[j])
            row.append(corr)
        corr_mtx.append(row)
    return corr_mtx


config = configs.PV_config
root_address = '/home/ssere004/Malaria/merge/extractor/' + config['organism_name']

meth_address_list = glob.glob(root_address + "/*CX_report.txt")
print(meth_address_list)
seq_address = '/home/ssere004/Malaria/genome_assembelies/Plasmodium_vivax/V1/PlasmoDB-48_PvivaxP01_Genome.fasta'

print(get_pearson_mtx(meth_address_list, seq_address, config))

