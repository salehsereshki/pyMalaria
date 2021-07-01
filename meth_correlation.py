import methods
import numpy as np
from scipy.stats import pearsonr
import input_parser as input_parser
import constants as constants

def get_pearson_mtx(meht_address_lst, seq_address):
    sample_meths = []
    sequences = input_parser.readfasta(seq_address)
    chros = list(sequences.keys())
    for i in meht_address_lst:
        methylations = input_parser.read_methylations(i)
        meth_seq , context_seq = input_parser.make_meth_string(methylations, sequences, 1)
        meth_str = []
        for chro in (chros):
            meth_str = meth_str + meth_seq[chro]
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
