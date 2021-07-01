from Bio import SeqIO
import pandas as pd
import numpy as np
import input_parser as input_parser
import constants as contants





def get_coverage(methylations):
    meth_sum = methylations['meth'].sum()
    unmeth_sum = methylations['unmeth'].sum()
    print((meth_sum + unmeth_sum)/len(methylations))



def exon_meth_counts(exon_df, meth_seq):
    exon_meth_p = []
    exon_meth_n = []
    flac_up_p = []
    flac_up_n = []
    flac_down_p = []
    flac_down_n = []
    for index, row in exon_df.iterrows():
        start = row['start']
        end = row['end']
        strand = row['strand']
        fp_start = end
        fp_end = end + 150
        fn_start = start - 150
        fn_end = start
        if strand == '-':
            fp_start = start - 150
            fp_end = start
            fn_start = end
            fn_end = end + 150
        meths_p = []
        meths_n = []
        for i in range(end - start):
            nucl_meth = float(meth_seq[row['chr']][i])
            if nucl_meth > 0:
                meths_p.append(nucl_meth)
                meths_n.append(0)
            if nucl_meth < 0:
                meths_n.append(-1*nucl_meth)
                meths_p.append(0)
            if nucl_meth == 0:
                meths_p.append(nucl_meth)
                meths_n.append(nucl_meth)
        meth_up_p = []
        meth_up_n = []
        for i in range(fp_end - fp_start):
            nucl_meth = float(meth_seq[row['chr']][i])
            if nucl_meth > 0:
                meth_up_p.append(nucl_meth)
                meth_up_n.append(0)
            if nucl_meth < 0:
                meth_up_n.append(-1*nucl_meth)
                meth_up_p.append(0)
            if nucl_meth == 0:
                meth_up_p.append(nucl_meth)
                meth_up_n.append(nucl_meth)
        meth_down_p = []
        meth_down_n = []
        for i in range(fn_end - fn_start):
            nucl_meth = float(meth_seq[row['chr']][i])
            if nucl_meth > 0:
                meth_down_p.append(nucl_meth)
                meth_down_n.append(0)
            if nucl_meth < 0:
                meth_down_n.append(-1*nucl_meth)
                meth_down_p.append(0)
            if nucl_meth == 0:
                meth_down_p.append(nucl_meth)
                meth_down_n.append(nucl_meth)
        exon_meth_p.append(meths_p)
        exon_meth_n.append(meths_n)
        flac_up_p.append(meth_up_p)
        flac_up_n.append(meth_up_n)
        flac_down_p.append(meth_down_p)
        flac_down_n.append(meth_down_n)
    return exon_meth_p, exon_meth_n, flac_up_p, flac_up_n, flac_down_p, flac_down_n

def make_exon_average_meths(exon_meth_p, exon_meth_n, flac_up_p, flac_up_n, flac_down_p, flac_down_n):
    exon_meth_start_p =[]
    exon_meth_start_n = []
    exon_meth_end_p = []
    exon_meth_end_n = []

    flac_up_avg_p = []
    flac_up_avg_n = []

    flac_down_avg_p = []
    flac_down_avg_n = []
    for i in range(150):
        sum_p = 0
        count_p = 0
        sum_n = 0
        count_n = 0
        for j in range(len(exon_meth_p)):
            if i < len(exon_meth_p[j]):
                sum_p += get_zero_instead(exon_meth_p[j][i])
                count_p += 1
            if i < len(exon_meth_n[j]):
                sum_n += get_zero_instead(exon_meth_n[j][i])
                count_n += 1
        exon_meth_start_p.append(float(sum_p)/ count_p)
        exon_meth_start_n.append(float(sum_n)/ count_n)
    for i in range(150):
        sum_p = 0
        count_p = 0
        sum_n = 0
        count_n = 0
        for j in range(len(exon_meth_p)):
            pstion_p = len(exon_meth_p[j]) - i - 1
            if pstion_p > 0:
                sum_p += get_zero_instead(exon_meth_p[j][pstion_p])
                count_p += 1
            pstion_n = len(exon_meth_p[j]) - i - 1
            if pstion_n > 0:
                sum_n += get_zero_instead(exon_meth_n[j][pstion_n])
                count_n += 1
        exon_meth_end_p.append(float(sum_p)/ count_p)
        exon_meth_end_n.append(float(sum_n)/ count_n)
    for i in range(150):
        sum_p = 0
        count_p = 0
        sum_n = 0
        count_n = 0
        for j in range(len(flac_up_p)):
            sum_p += get_zero_instead(flac_up_p[j][i])
            count_p += 1
            sum_n += get_zero_instead(flac_up_n[j][i])
            count_n += 1
        flac_up_avg_p.append(float(sum_p)/ count_p)
        flac_up_avg_n.append(float(sum_n)/ count_n)
    for i in range(150):
        sum_p = 0
        count_p = 0
        sum_n = 0
        count_n = 0
        for j in range(len(flac_down_p)):
            sum_p += get_zero_instead(flac_down_p[j][i])
            count_p += 1
            sum_n += get_zero_instead(flac_down_n[j][i])
            count_n += 1
        flac_down_avg_p.append(float(sum_p) / count_p)
        flac_down_avg_n.append(float(sum_n) / count_n)

    return exon_meth_start_p, exon_meth_start_n, exon_meth_end_p, exon_meth_end_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n


def get_zero_instead(mth):
    if mth == contants.NON_METH_TAG:
        return 0
    else:
        return mth

def get_exon_output(exon_meth_start_p, exon_meth_start_n, exon_meth_end_p, exon_meth_end_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n):
    down_p = flac_down_avg_p + exon_meth_start_p
    down_n = flac_down_avg_n + exon_meth_start_n

    up_p = exon_meth_end_p + flac_up_avg_p
    up_n = exon_meth_end_n + flac_up_avg_n

    return down_p, down_n, up_p, up_n




def get_counts(arr):
    counts_p = np.zeros(300)
    counts_n = np.zeros(300)

    for count_indx in range(len(counts_n)):
        for j in range(len(arr)):
            if float(arr[j][0][count_indx]) > 0:
                counts_p[count_indx] += 1
            if float(arr[j][0][count_indx]) < 0:
                counts_n[count_indx] += 1
    return counts_p[151], counts_n[151]







