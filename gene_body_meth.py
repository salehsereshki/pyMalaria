import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import pyMalaria.gbm_computations as gbmc

def get_gene_meth(meth_seq, genes_df,  bin_num, threshold = 0.1, flanking_size = 2000):
    genes_avg_p = np.zeros(bin_num, dtype=np.double)
    genes_avg_n = np.zeros(bin_num, dtype=np.double)
    flac_up_avg_p = np.zeros(bin_num, dtype=np.double)
    flac_up_avg_n = np.zeros(bin_num, dtype=np.double)
    flac_down_avg_p = np.zeros(bin_num, dtype=np.double)
    flac_down_avg_n = np.zeros(bin_num, dtype=np.double)

    gene_bins_sum = np.zeros(bin_num)
    flac_up_bins_sum = np.zeros(bin_num)
    flac_down_bins_sum = np.zeros(bin_num)

    for index, row in genes_df.iterrows():
        seq = meth_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
        #prev_gene_end = 0
        #next_gene_start = sys.maxsize
        #if index != 0 and (genes_df.iloc[index-1]['chr'] == row['chr']):
        #    prev_gene_end = genes_df[index - 1]['end']
        #if index != len(genes_df) and genes_df.iloc[index+1]['chr'] == row['chr']:
        #    next_gene_start = genes_df[index + 1]['start']

        check_flac_down_overlap = (index == 0 or (genes_df.iloc[index-1]['chr'] == row['chr'] and genes_df.iloc[index - 1]['end'] < row['start'] - flanking_size))
        check_flac_up_overlap = (index == len(genes_df) - 1 or (genes_df.iloc[index+1]['chr'] == row['chr'] and genes_df.iloc[index + 1]['start'] > row['end'] + flanking_size))

        #if row['strand'] == '-':
        #    seq = seq[::-1]
        #    check_flac_down_overlap, check_flac_up_overlap = check_flac_up_overlap, check_flac_down_overlap
            #prev_gene_end, next_gene_start = next_gene_start, prev_gene_end

        flac_bin_size = int(flanking_size/bin_num)
        gene_bin_size = int((row['end'] - row['start'])/bin_num)
        for i in range(bin_num):
            m_p, m_n = get_meth_percentage(seq[i*flac_bin_size: (i+1) * flac_bin_size], threshold)
            if m_n != None and m_p != None and check_flac_down_overlap:
                if index == 0 or (genes_df.iloc[index-1]['chr'] == row['chr'] and genes_df.iloc[index - 1]['end'] < row['start'] - flanking_size):
                    flac_down_avg_p[i] = update_mean(flac_down_avg_p[i], flac_down_bins_sum[i], m_p, flac_bin_size)
                    flac_down_avg_n[i] = update_mean(flac_down_avg_n[i], flac_down_bins_sum[i], m_n, flac_bin_size)
                    flac_down_bins_sum[i] += flac_bin_size

            m_p, m_n = get_meth_percentage(seq[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size], threshold)
            if m_n != None and m_p != None:
                genes_avg_p[i] = update_mean(genes_avg_p[i], gene_bins_sum[i], m_p, gene_bin_size)
                genes_avg_n[i] = update_mean(genes_avg_n[i], gene_bins_sum[i], m_n, gene_bin_size)
                gene_bins_sum[i] += gene_bin_size

            m_p, m_n = get_meth_percentage(seq[i*flac_bin_size + len(seq) - flanking_size: (i+1) * flac_bin_size + len(seq) - flanking_size], threshold)
            if m_n != None and m_p != None and check_flac_up_overlap:
                flac_up_avg_p[i] = update_mean(flac_up_avg_p[i], flac_up_bins_sum[i], m_p, flac_bin_size)
                flac_up_avg_n[i] = update_mean(flac_up_avg_n[i], flac_up_bins_sum[i], m_n, flac_bin_size)
                flac_up_bins_sum[i] += flac_bin_size

    return genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n


def get_gene_meth_count_based(meth_count_seq, unmeth_count_seq,  genes_df,  bin_num, flanking_size = 2000):
    genes_avg_p = np.zeros(bin_num, dtype=np.double)
    genes_avg_n = np.zeros(bin_num, dtype=np.double)
    flac_up_avg_p = np.zeros(bin_num, dtype=np.double)
    flac_up_avg_n = np.zeros(bin_num, dtype=np.double)
    flac_down_avg_p = np.zeros(bin_num, dtype=np.double)
    flac_down_avg_n = np.zeros(bin_num, dtype=np.double)

    gene_bins_sum = np.zeros(bin_num)
    flac_up_bins_sum = np.zeros(bin_num)
    flac_down_bins_sum = np.zeros(bin_num)

    count_flac_down_overlap = 0
    count_flac_up_overlap = 0
    count_all = 0

    for index, row in genes_df.iterrows():
        prev_gene_end = -1
        next_gene_start = sys.maxsize

        if index > 0 and genes_df.iloc[index-1]['chr'] == row['chr']:
            prev_gene_end = genes_df.iloc[index-1]['end']
        if index < len(genes_df)-1 and genes_df.iloc[index+1]['chr'] == row['chr']:
            next_gene_start = genes_df.iloc[index+1]['start']

        if False:
            start_seq_position = max(row['start'] - flanking_size, prev_gene_end)
            end_seq_position = min(row['end'] + flanking_size, next_gene_start)

            seq_meth = meth_count_seq[row['chr']][start_seq_position: end_seq_position]
            seq_unmeth = unmeth_count_seq[row['chr']][start_seq_position: end_seq_position]

            # if row['strand'] == '-':
            #     seq_meth = seq_meth[::-1]
            #     seq_unmeth = seq_unmeth[::-1]

            flac_bin_size = int(flanking_size/bin_num)
            gene_bin_size = int((row['end'] - row['start'])/bin_num)

            bin_start = row['start'] - start_seq_position - flac_bin_size
            a = bin_start
            i = bin_num - 1
            while i >= 0 and bin_start > prev_gene_end:
                m_p, m_n = get_meth_percentage_count_based(seq_meth[bin_start: bin_start + flac_bin_size], seq_unmeth[bin_start: bin_start + flac_bin_size])
                if m_n != None and m_p != None:
                    flac_down_avg_p[i] = update_mean(flac_down_avg_p[i], flac_down_bins_sum[i], m_p, flac_bin_size)
                    flac_down_avg_n[i] = update_mean(flac_down_avg_n[i], flac_down_bins_sum[i], m_n, flac_bin_size)
                    flac_down_bins_sum[i] += flac_bin_size
                i = i - 1
                bin_start = bin_start - flac_bin_size
                if count_all == 0:
                    print(bin_start)

            bin_start = row['start'] - start_seq_position
            b = bin_start
            i = 0
            while i < bin_num:
                m_p, m_n = get_meth_percentage_count_based(seq_meth[bin_start: bin_start + gene_bin_size], seq_unmeth[bin_start: bin_start + gene_bin_size])
                if m_n != None and m_p != None:
                    genes_avg_p[i] = update_mean(genes_avg_p[i], gene_bins_sum[i], m_p, gene_bin_size)
                    genes_avg_n[i] = update_mean(genes_avg_n[i], gene_bins_sum[i], m_n, gene_bin_size)
                    gene_bins_sum[i] += gene_bin_size
                i = i + 1
                bin_start = bin_start + gene_bin_size
                if count_all == 0:
                    print(bin_start)

            bin_start = len(seq_meth) - (end_seq_position - row['end'])
            c = bin_start
            i = 0
            while i < bin_num and bin_start + flac_bin_size < next_gene_start:
                m_p, m_n = get_meth_percentage_count_based(seq_meth[bin_start: bin_start + flac_bin_size], seq_unmeth[bin_start: bin_start + flac_bin_size])
                if m_n != None and m_p != None:
                    flac_up_avg_p[i] = update_mean(flac_up_avg_p[i], flac_up_bins_sum[i], m_p, flac_bin_size)
                    flac_up_avg_n[i] = update_mean(flac_up_avg_n[i], flac_up_bins_sum[i], m_n, flac_bin_size)
                    flac_up_bins_sum[i] += flac_bin_size
                i = i + 1
                bin_start = bin_start + flac_bin_size
                if count_all == 0:
                    print(bin_start)
            if count_all < 1:
                print(a, b, c, row['start'], row['end'], prev_gene_end, next_gene_start, len(seq_meth), len(seq_unmeth))
                count_all+=1


        else:
            seq_meth = meth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
            seq_unmeth = unmeth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
            flac_bin_size = int(flanking_size/bin_num)
            gene_bin_size = int((row['end'] - row['start'])/bin_num)
            flac_down_num_binsoverlap = min((row['start'] - prev_gene_end)/flac_bin_size, bin_num)
            flac_up_num_binsoverlap = min((next_gene_start - row['end'])/flac_bin_size, bin_num)

            check_flac_down_overlap = (index == 0 or (genes_df.iloc[index-1]['chr'] == row['chr'] and genes_df.iloc[index - 1]['end'] < row['start'] - flanking_size))
            check_flac_up_overlap = (index == len(genes_df) - 1 or (genes_df.iloc[index+1]['chr'] == row['chr'] and genes_df.iloc[index + 1]['start'] > row['end'] + flanking_size))

            if check_flac_down_overlap:
                count_flac_down_overlap += 1
            if check_flac_up_overlap:
                count_flac_up_overlap += 1
            count_all += 1

            if row['strand'] == '-':
                seq_meth = seq_meth[::-1]
                seq_unmeth = seq_unmeth[::-1]
                check_flac_down_overlap, check_flac_up_overlap = check_flac_up_overlap, check_flac_down_overlap
                flac_down_num_binsoverlap, flac_up_num_binsoverlap = flac_up_num_binsoverlap, flac_down_num_binsoverlap

            for i in range(bin_num):
                m_p, m_n = get_meth_percentage_count_based(seq_meth[i*flac_bin_size: (i+1) * flac_bin_size], seq_unmeth[i*flac_bin_size: (i+1) * flac_bin_size])
                if m_n != None and m_p != None and i > bin_num - flac_down_num_binsoverlap - 1:
                    flac_down_avg_p[i] = update_mean(flac_down_avg_p[i], flac_down_bins_sum[i], m_p, flac_bin_size)
                    flac_down_avg_n[i] = update_mean(flac_down_avg_n[i], flac_down_bins_sum[i], m_n, flac_bin_size)
                    flac_down_bins_sum[i] += flac_bin_size

                m_p, m_n = get_meth_percentage_count_based(seq_meth[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size], seq_unmeth[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size])
                if m_n != None and m_p != None:
                    genes_avg_p[i] = update_mean(genes_avg_p[i], gene_bins_sum[i], m_p, gene_bin_size)
                    genes_avg_n[i] = update_mean(genes_avg_n[i], gene_bins_sum[i], m_n, gene_bin_size)
                    gene_bins_sum[i] += gene_bin_size

                m_p, m_n = get_meth_percentage_count_based(seq_meth[i*flac_bin_size + len(seq_meth) - flanking_size: (i+1) * flac_bin_size + len(seq_meth) - flanking_size], seq_unmeth[i*flac_bin_size + len(seq_unmeth) - flanking_size: (i+1) * flac_bin_size + len(seq_unmeth) - flanking_size])
                if m_n != None and m_p != None and i < flac_up_num_binsoverlap:
                    flac_up_avg_p[i] = update_mean(flac_up_avg_p[i], flac_up_bins_sum[i], m_p, flac_bin_size)
                    flac_up_avg_n[i] = update_mean(flac_up_avg_n[i], flac_up_bins_sum[i], m_n, flac_bin_size)
                    flac_up_bins_sum[i] += flac_bin_size
    print(count_flac_down_overlap, count_flac_up_overlap, count_all)
    return genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n


def update_mean(mean, sum_weights, new_value, weight):
    return ((mean * sum_weights) + (new_value * weight)) / (sum_weights + weight)




def get_gene_meth_count_based_average_bin(meth_count_seq, unmeth_count_seq,  genes_df,  bin_num, flanking_size = 2000):
    genes_sumcounts_p = np.zeros([2, bin_num], dtype=np.double)
    genes_sumcounts_n = np.zeros([2, bin_num], dtype=np.double)
    flac_up_sumcounts_p = np.zeros([2, bin_num], dtype=np.double)
    flac_up_sumcounts_n = np.zeros([2, bin_num], dtype=np.double)
    flac_down_sumcounts_p = np.zeros([2, bin_num], dtype=np.double)
    flac_down_sumcounts_n = np.zeros([2, bin_num], dtype=np.double)


    for index, row in genes_df.iterrows():
        seq_meth = meth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
        seq_unmeth = unmeth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
        check_flac_down_overlap = (index == 0 or (genes_df.iloc[index-1]['chr'] == row['chr'] and genes_df.iloc[index - 1]['end'] < row['start'] - flanking_size))
        check_flac_up_overlap = (index == len(genes_df) - 1 or (genes_df.iloc[index+1]['chr'] == row['chr'] and genes_df.iloc[index + 1]['start'] > row['end'] + flanking_size))

        #if row['strand'] == '-':
            #seq_meth = seq_meth[::-1]
            #seq_unmeth = seq_unmeth[::-1]
            #check_flac_down_overlap, check_flac_up_overlap = check_flac_up_overlap, check_flac_down_overlap
        flac_bin_size = int(flanking_size/bin_num)
        gene_bin_size = int((row['end'] - row['start'])/bin_num)
        for i in range(bin_num):
            if check_flac_down_overlap:
                sm = seq_meth[i*flac_bin_size: (i+1) * flac_bin_size]
                sum = seq_unmeth[i*flac_bin_size: (i+1) * flac_bin_size]

                flac_down_sumcounts_p[0][i] += np.sum(sm[sm > 0])
                flac_down_sumcounts_p[1][i] += np.sum(sum[sum > 0])

                flac_down_sumcounts_n[0][i] += np.sum(sm[sm < 0])
                flac_down_sumcounts_n[1][i] += np.sum(sum[sum < 0])

            sm = seq_meth[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size]
            sum = seq_unmeth[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size]

            genes_sumcounts_p[0][i] += np.sum(sm[sm > 0])
            genes_sumcounts_p[1][i] += np.sum(sum[sum > 0])

            genes_sumcounts_n[0][i] += np.sum(sm[sm < 0])
            genes_sumcounts_n[1][i] += np.sum(sum[sum < 0])

            if check_flac_up_overlap:
                sm = seq_meth[i*flac_bin_size + len(seq_meth) - flanking_size: (i+1) * flac_bin_size + len(seq_meth) - flanking_size]
                sum = seq_unmeth[i*flac_bin_size + len(seq_unmeth) - flanking_size: (i+1) * flac_bin_size + len(seq_unmeth) - flanking_size]

                flac_up_sumcounts_p[0][i] += np.sum(sm[sm > 0])
                flac_up_sumcounts_p[1][i] += np.sum(sum[sum > 0])

                flac_up_sumcounts_n[0][i] += np.sum(sm[sm < 0])
                flac_up_sumcounts_n[1][i] += np.sum(sum[sum < 0])


    return avg_count(genes_sumcounts_p, bin_num), avg_count(genes_sumcounts_n, bin_num),\
           avg_count(flac_up_sumcounts_p, bin_num), avg_count(flac_up_sumcounts_n, bin_num),\
           avg_count(flac_down_sumcounts_p, bin_num), avg_count(flac_down_sumcounts_n, bin_num)

def avg_count(sum_count, bin_num):
    avgs = np.zeros(bin_num, dtype=np.double)
    for i in range(bin_num):
        avgs[i] = sum_count[0][i] / (sum_count[0][i] + sum_count[1][i])
    return avgs

def strand_specific_meth(organism_name, methylations, genes_seq, threshold, coverage_threshold):
    genes_same_C = 0
    genes_same_meC = 0
    genes_non_same_C = 0
    genes_non_same_meC = 0

    nongenes_p_C = 0
    nongenes_p_meC = 0
    nongenes_n_C = 0
    nongenes_n_meC = 0

    for index, row in methylations.iterrows():
        gene_status = genes_seq[row['chr']][row['position'] - 1]
        if row['meth']+row['unmeth'] > coverage_threshold:
            is_meth = row['meth']/(row['meth']+row['unmeth']) > threshold
        else:
            continue
        c_strand = row['strand']

        if gene_status != 0:
            if (c_strand == '+' and gene_status == 1) or (c_strand == '-' and gene_status == -1):
                genes_same_C+=1
                if is_meth:
                    genes_same_meC+=1
            else:
                genes_non_same_C+=1
                if is_meth:
                    genes_non_same_meC+=1
        else:
            if c_strand == '+':
                nongenes_p_C+=1
                if is_meth:
                    nongenes_p_meC+=1
            else:
                nongenes_n_C+=1
                if is_meth:
                    nongenes_n_meC+=1

    res = [genes_same_meC/genes_same_C, genes_non_same_meC/genes_non_same_C, nongenes_p_meC/nongenes_p_C, nongenes_n_meC/nongenes_n_C]

    xloc = range(4)
    labels = ['gene sense\n strand', 'gene anti sense \n strand', 'non gene \n positive strand', 'non gene \n negative strand']
    barWidth = .5

    plt.bar(xloc, res, width=barWidth)

    plt.xticks(xloc, labels)
    plt.savefig('gene_nongene_strand' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()

    return


def get_gene_body_meth_half(meth_seq, genes_df, organism_name):

    first_half_meth_p_Cp = np.full(len(genes_df), np.NaN)
    first_half_meth_n_Cp = np.full(len(genes_df), np.NaN)

    first_half_meth_p_Cn = np.full(len(genes_df), np.NaN)
    first_half_meth_n_Cn = np.full(len(genes_df), np.NaN)

    second_half_meth_p_Cp = np.full(len(genes_df), np.NaN)
    second_half_meth_n_Cp = np.full(len(genes_df), np.NaN)

    second_half_meth_p_Cn = np.full(len(genes_df), np.NaN)
    second_half_meth_n_Cn = np.full(len(genes_df), np.NaN)

    res_input = [first_half_meth_p_Cp,first_half_meth_n_Cp,
     first_half_meth_p_Cn, first_half_meth_n_Cn,
     second_half_meth_p_Cp, second_half_meth_n_Cp,
     second_half_meth_p_Cn, second_half_meth_n_Cn
     ]

    bins = np.full(len(genes_df), np.NaN)

    for index, row in genes_df.iterrows():
        middle = int(int(row['end'] - row['start'])/2 + row['start'])
        bins[index] = int(int(row['end'] - row['start'])/2)
        if row['strand'] == '-':
            first_half_meth_n_Cp[index], first_half_meth_n_Cn[index] = get_meth_percentage(meth_seq[row['chr']][row['start']: middle], 0.1)
            second_half_meth_n_Cp[index], second_half_meth_n_Cn[index] = get_meth_percentage(meth_seq[row['chr']][middle: row['end']], 0.1)
        else:
            first_half_meth_p_Cp[index], first_half_meth_p_Cn[index] = get_meth_percentage(meth_seq[row['chr']][row['start']: middle], 0.1)
            second_half_meth_p_Cp[index], second_half_meth_p_Cn[index] = get_meth_percentage(meth_seq[row['chr']][middle: row['end']], 0.1)

    res_out = []
    for i in range(len(res_input)):
        res_out.append(np.nansum(np.multiply(bins, res_input[i])) / np.sum(bins[np.logical_not(np.isnan(res_input[i]))]))

    res_out_Cp = [res_out[0], res_out[1], res_out[4], res_out[5]]
    res_out_Cn = [res_out[2], res_out[3], res_out[6], res_out[7]]

    xloc = range(4)
    labels = ['1st g+', '1st g-', '2nd g+', '2nd g-']
    barWidth = .5

    p1 = plt.bar(xloc, res_out_Cp, width=barWidth)
    p2 = plt.bar(xloc, res_out_Cn, bottom=res_out_Cp, width=barWidth)

    plt.xticks(xloc, labels)
    plt.legend((p1[0], p2[0]), ('C +', 'C -'))
    plt.savefig('genebody_strand_half' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()


    return res_out





def get_meth_percentage(meth_stat, threshold):
    countCs_p = 0
    countCs_n = 0
    countMethCs_p = 0
    countMethCs_n = 0
    for i in range(len(meth_stat)):
        if float(meth_stat[i]) > 0:
            countCs_p += 1
        elif float(meth_stat[i]) < 0:
            countCs_n += 1
        if float(meth_stat[i]) > threshold:
            countMethCs_p += 1
        elif float(meth_stat[i]) < -1 * threshold:
            countMethCs_n += 1
    if countCs_p != 0 and countCs_n != 0:
        return float(countMethCs_p) / countCs_p, float(countMethCs_n) / countCs_n # QUESTION
    else:
        return None, None

def get_meth_percentage_count_based(meth_stat, unmeth_state):

    counts_meth_p = 0
    counts_unmeth_p = 0

    counts_meth_n = 0
    counts_unmeth_n = 0

    for i in meth_stat:
        if i > 0:
            counts_meth_p += i
        else:
            counts_meth_n += -1 * i

    for i in unmeth_state:
        if i > 0:
            counts_unmeth_p += i
        else:
            counts_unmeth_n += -1 * i

    if counts_meth_p+counts_unmeth_p != 0 and counts_meth_n+counts_unmeth_n != 0:
        return float(counts_meth_p)/ (counts_meth_p+counts_unmeth_p), float(counts_meth_n)/ (counts_meth_n+counts_unmeth_n)
    else:
        return None, None

def get_average(meth_profs_df, is_flanking, bin_num, sum_bin_sizes = None):

    meth_df = pd.DataFrame(meth_profs_df, columns = range(bin_num))
    res = []
    if not is_flanking:
        for i in range(len(meth_df.columns)):
            res.append(meth_df[meth_df.columns[i]].sum(skipna=True)/sum_bin_sizes[i])
    else:
        for i in range(len(meth_df.columns)):
            res.append(meth_df[meth_df.columns[i]].mean(skipna=True))
    return res



def plot_gene_body_meth(organism_name, meth_seq, genes_df, bin_num, threshold = 0.1):
    if len(meth_seq) == 2:
        genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n = get_gene_meth_count_based(meth_seq[0], meth_seq[1], genes_df,  bin_num)
    else:
        genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n = get_gene_meth(meth_seq, genes_df,  bin_num, threshold = threshold)

    final_p = np.concatenate((flac_down_avg_p , genes_avg_p , flac_up_avg_p))
    final_n = np.concatenate((flac_down_avg_n , genes_avg_n , flac_up_avg_n))
    print(final_p)
    print(final_n)
    # yticks = [0]
    # ylabels = ['']
    # plt.yticks(yticks, ylabels)
    #plt.tick_params(left=False, labelleft=False)
    plt.box(False)
    plt.ylabel("$meC/C$")
    ticks = [0, bin_num, bin_num * 2, bin_num * 3]
    labels = ['      5\' flanking region', '           gene body', '       3\' flanking region', '']
    plt.xticks(ticks, labels, horizontalalignment='left')
    # plt.tick_params(axis='x', colors='black', direction='out', length=10, width=1,  pad = 4)
    plt.grid(False)
    plt.style.use('seaborn')
    plt.plot(range(0, 3 * bin_num), final_p, color='blue', linewidth=4.0)
    plt.plot(range(0, 3 * bin_num), final_n, color='red', linewidth=4.0)
    plt.axhline(y=0.0, color='black', linestyle='-')
    plt.axvline(x=0.0, color='black', linestyle='-')
    plt.rcParams['axes.facecolor'] = 'white'
    plt.savefig('genebody_' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()

def plot_template_non_template(organism_name, meth_seq, genes_df, threshold=0.1):
    temp_meth, nontemp_meth = gbmc.template_methylation(meth_seq, genes_df, threshold=threshold)
    res = [temp_meth, nontemp_meth]
    xloc = range(2)
    labels = ['temp', 'non-temp']
    barWidth = .2

    plt.bar(xloc, res, width=barWidth)

    plt.xticks(xloc, labels)
    plt.savefig('temp-nontemp-' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()
