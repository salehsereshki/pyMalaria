import numpy as np
import sys

#Count_based, average_bin, method1
def get_gene_meth_count_based_average_bin(meth_count_seq, unmeth_count_seq,  genes_df,  bin_num, flanking_size = 2000):
    print('1')
    genes_sumcounts_p = np.zeros([2, bin_num], dtype=np.double)
    genes_sumcounts_n = np.zeros([2, bin_num], dtype=np.double)
    flac_up_sumcounts_p = np.zeros([2, bin_num], dtype=np.double)
    flac_up_sumcounts_n = np.zeros([2, bin_num], dtype=np.double)
    flac_down_sumcounts_p = np.zeros([2, bin_num], dtype=np.double)
    flac_down_sumcounts_n = np.zeros([2, bin_num], dtype=np.double)


    for index, row in genes_df.iterrows():
        prev_gene_end = -1
        next_gene_start = sys.maxsize

        if index > 0 and genes_df.iloc[index-1]['chr'] == row['chr']:
            prev_gene_end = genes_df.iloc[index-1]['end']
        if index < len(genes_df)-1 and genes_df.iloc[index+1]['chr'] == row['chr']:
            next_gene_start = genes_df.iloc[index+1]['start']


        start_seq_position = max(row['start'] - flanking_size, prev_gene_end)
        end_seq_position = min(row['end'] + flanking_size, next_gene_start)

        seq_meth = meth_count_seq[row['chr']][start_seq_position: end_seq_position]
        seq_unmeth = unmeth_count_seq[row['chr']][start_seq_position: end_seq_position]

        if row['strand'] == '-':
            seq_meth = seq_meth[::-1]
            seq_unmeth = seq_unmeth[::-1]

        flac_bin_size = int(flanking_size/bin_num)
        gene_bin_size = int((row['end'] - row['start'])/bin_num)

        if row['strand'] == '+':
            bin_start = row['start'] - start_seq_position - flac_bin_size
        else:
            bin_start = end_seq_position - row['end'] - flac_bin_size

        i = bin_num - 1
        while i >= 0 and bin_start >= 0:
            sm = seq_meth[bin_start: bin_start + flac_bin_size]
            sum = seq_unmeth[bin_start: bin_start + flac_bin_size]
            flac_down_sumcounts_p[0][i] += np.sum(sm[sm > 0])
            flac_down_sumcounts_p[1][i] += np.sum(sum[sum > 0])
            flac_down_sumcounts_n[0][i] += abs(np.sum(sm[sm < 0]))
            flac_down_sumcounts_n[1][i] += abs(np.sum(sum[sum < 0]))
            i = i - 1
            bin_start = bin_start - flac_bin_size

        if row['strand'] == '+':
            bin_start = row['start'] - start_seq_position
        else:
            bin_start = end_seq_position - row['end']
        i = 0
        while i < bin_num:
            sm = seq_meth[bin_start: bin_start + gene_bin_size]
            sum = seq_unmeth[bin_start: bin_start + gene_bin_size]
            genes_sumcounts_p[0][i] += np.sum(sm[sm > 0])
            genes_sumcounts_p[1][i] += np.sum(sum[sum > 0])
            genes_sumcounts_n[0][i] += abs(np.sum(sm[sm < 0]))
            genes_sumcounts_n[1][i] += abs(np.sum(sum[sum < 0]))
            i = i + 1
            bin_start = bin_start + gene_bin_size

        if row['strand'] == '+':
            bin_start = len(seq_meth) - (end_seq_position - row['end'])
        else:
            bin_start = len(seq_meth) - (row['start'] - start_seq_position)
        i = 0
        while i < bin_num and bin_start + flac_bin_size <= len(seq_meth):
            sm = seq_meth[bin_start: bin_start + flac_bin_size]
            sum = seq_unmeth[bin_start: bin_start + flac_bin_size]

            flac_up_sumcounts_p[0][i] += np.sum(sm[sm > 0])
            flac_up_sumcounts_p[1][i] += np.sum(sum[sum > 0])

            flac_up_sumcounts_n[0][i] += abs(np.sum(sm[sm < 0]))
            flac_up_sumcounts_n[1][i] += abs(np.sum(sum[sum < 0]))
            i = i + 1
            bin_start = bin_start + flac_bin_size



    return avg_count(genes_sumcounts_p, bin_num), avg_count(genes_sumcounts_n, bin_num),\
           avg_count(flac_up_sumcounts_p, bin_num), avg_count(flac_up_sumcounts_n, bin_num),\
           avg_count(flac_down_sumcounts_p, bin_num), avg_count(flac_down_sumcounts_n, bin_num)

#Count_base, flac_num_bin_overlap, method2
def get_gene_meth_count_based1(meth_count_seq, unmeth_count_seq,  genes_df,  bin_num, flanking_size = 2000):
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
        prev_gene_end = -1
        next_gene_start = sys.maxsize

        if index > 0 and genes_df.iloc[index-1]['chr'] == row['chr']:
            prev_gene_end = genes_df.iloc[index-1]['end']
        if index < len(genes_df)-1 and genes_df.iloc[index+1]['chr'] == row['chr']:
            next_gene_start = genes_df.iloc[index+1]['start']

        seq_meth = meth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
        seq_unmeth = unmeth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
        flac_bin_size = int(flanking_size/bin_num)
        gene_bin_size = int((row['end'] - row['start'])/bin_num)
        flac_down_num_binsoverlap = min((row['start'] - prev_gene_end)/flac_bin_size, bin_num)
        flac_up_num_binsoverlap = min((next_gene_start - row['end'])/flac_bin_size, bin_num)


        if row['strand'] == '-':
            seq_meth = seq_meth[::-1]
            seq_unmeth = seq_unmeth[::-1]
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
    return genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n

#Count_base, Start_bin, Method 3
def get_gene_meth_count_based2(meth_count_seq, unmeth_count_seq,  genes_df,  bin_num, flanking_size = 2000):
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

        test_flac_down = np.ones(5)
        test_gene = np.ones(5)
        test_flac_up = np.ones(5)


        if index > 0 and genes_df.iloc[index-1]['chr'] == row['chr']:
            prev_gene_end = genes_df.iloc[index-1]['end']
        if index < len(genes_df)-1 and genes_df.iloc[index+1]['chr'] == row['chr']:
            next_gene_start = genes_df.iloc[index+1]['start']

        start_seq_position = max(row['start'] - flanking_size, prev_gene_end)
        end_seq_position = min(row['end'] + flanking_size, next_gene_start)

        seq_meth = meth_count_seq[row['chr']][start_seq_position: end_seq_position]
        seq_unmeth = unmeth_count_seq[row['chr']][start_seq_position: end_seq_position]

        if row['strand'] == '-':
            seq_meth = seq_meth[::-1]
            seq_unmeth = seq_unmeth[::-1]

        flac_bin_size = int(flanking_size/bin_num)
        gene_bin_size = int((row['end'] - row['start'])/bin_num)

        if row['strand'] == '+':
            bin_start = row['start'] - start_seq_position - flac_bin_size
        else:
            bin_start = end_seq_position - row['end'] - flac_bin_size
        i = bin_num - 1
        while i >= 0 and bin_start >= 0:
            m_p, m_n = get_meth_percentage_count_based(seq_meth[bin_start: bin_start + flac_bin_size], seq_unmeth[bin_start: bin_start + flac_bin_size])
            if m_n != None and m_p != None:
                flac_down_avg_p[i] = update_mean(flac_down_avg_p[i], flac_down_bins_sum[i], m_p, flac_bin_size)
                flac_down_avg_n[i] = update_mean(flac_down_avg_n[i], flac_down_bins_sum[i], m_n, flac_bin_size)
                flac_down_bins_sum[i] += flac_bin_size
            test_flac_down[i] = bin_start
            i = i - 1
            bin_start = bin_start - flac_bin_size

        if row['strand'] == '+':
            bin_start = row['start'] - start_seq_position
        else:
            bin_start = end_seq_position - row['end']
        i = 0
        while i < bin_num:
            m_p, m_n = get_meth_percentage_count_based(seq_meth[bin_start: bin_start + gene_bin_size], seq_unmeth[bin_start: bin_start + gene_bin_size])
            if m_n != None and m_p != None:
                genes_avg_p[i] = update_mean(genes_avg_p[i], gene_bins_sum[i], m_p, gene_bin_size)
                genes_avg_n[i] = update_mean(genes_avg_n[i], gene_bins_sum[i], m_n, gene_bin_size)
                gene_bins_sum[i] += gene_bin_size
            test_gene[i] = bin_start
            i = i + 1
            bin_start = bin_start + gene_bin_size


        if row['strand'] == '+':
            bin_start = len(seq_meth) - (end_seq_position - row['end'])
        else:
            bin_start = len(seq_meth) - (row['start'] - start_seq_position)
        i = 0
        while i < bin_num and bin_start + flac_bin_size <= len(seq_meth):
            m_p, m_n = get_meth_percentage_count_based(seq_meth[bin_start: bin_start + flac_bin_size], seq_unmeth[bin_start: bin_start + flac_bin_size])
            if m_n != None and m_p != None:
                flac_up_avg_p[i] = update_mean(flac_up_avg_p[i], flac_up_bins_sum[i], m_p, flac_bin_size)
                flac_up_avg_n[i] = update_mean(flac_up_avg_n[i], flac_up_bins_sum[i], m_n, flac_bin_size)
                flac_up_bins_sum[i] += flac_bin_size
            test_flac_up[i] = bin_start
            i = i + 1
            bin_start = bin_start + flac_bin_size

    print(count_flac_down_overlap, count_flac_up_overlap, count_all)
    return genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n

#CLevel, flac_num_bin_overlap, method4
def get_gene_meth(meth_seq, genes_df,  bin_num, exp_prcnt, threshold = 0.1, flanking_size = 500):
    genes_avg_p = np.zeros(bin_num, dtype=np.double)
    genes_avg_n = np.zeros(bin_num, dtype=np.double)
    flac_up_avg_p = np.zeros(bin_num, dtype=np.double)
    flac_up_avg_n = np.zeros(bin_num, dtype=np.double)
    flac_down_avg_p = np.zeros(bin_num, dtype=np.double)
    flac_down_avg_n = np.zeros(bin_num, dtype=np.double)

    gene_bins_sum = np.zeros(bin_num)
    flac_up_bins_sum = np.zeros(bin_num)
    flac_down_bins_sum = np.zeros(bin_num)

    me_df_t = np.zeros((len(genes_df), 3 * bin_num))
    me_df_nt = np.zeros((len(genes_df), 3 * bin_num))


    for index, row in genes_df.iterrows():
        is_template = row['strand'] == '+'
        prev_gene_end = -1
        next_gene_start = sys.maxsize

        if index > 0 and genes_df.iloc[index-1]['chr'] == row['chr']:
            prev_gene_end = genes_df.iloc[index-1]['end']
        if index < len(genes_df)-1 and genes_df.iloc[index+1]['chr'] == row['chr']:
            next_gene_start = genes_df.iloc[index+1]['start']

        seq_meth = meth_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]

        flac_bin_size = int(flanking_size/bin_num)
        gene_bin_size = int((row['end'] - row['start'])/bin_num)
        flac_down_num_binsoverlap = min((row['start'] - prev_gene_end)/flac_bin_size, bin_num)
        flac_up_num_binsoverlap = min((next_gene_start - row['end'])/flac_bin_size, bin_num)

        if not is_template:
            seq_meth = seq_meth[::-1]
            flac_down_num_binsoverlap, flac_up_num_binsoverlap = flac_up_num_binsoverlap, flac_down_num_binsoverlap

        for i in range(bin_num):
            m_p, m_n = get_meth_percentage(seq_meth[i*flac_bin_size: (i+1) * flac_bin_size], threshold)
            if not is_template:
                m_p, m_n = m_n, m_p
            me_df_t[index][i] = m_p
            me_df_nt[index][i] = m_n
            if m_n != None and m_p != None and i > bin_num - flac_down_num_binsoverlap - 1:
                flac_down_avg_p[i] = update_mean(flac_down_avg_p[i], flac_down_bins_sum[i], m_p, flac_bin_size)
                flac_down_avg_n[i] = update_mean(flac_down_avg_n[i], flac_down_bins_sum[i], m_n, flac_bin_size)
                flac_down_bins_sum[i] += flac_bin_size

            m_p, m_n = get_meth_percentage(seq_meth[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size], threshold)
            if not is_template:
                m_p, m_n = m_n, m_p
            me_df_t[index][i+5] = m_p
            me_df_nt[index][i+5] = m_n
            if m_n != None and m_p != None:
                genes_avg_p[i] = update_mean(genes_avg_p[i], gene_bins_sum[i], m_p, gene_bin_size)
                genes_avg_n[i] = update_mean(genes_avg_n[i], gene_bins_sum[i], m_n, gene_bin_size)
                gene_bins_sum[i] += gene_bin_size

            m_p, m_n = get_meth_percentage(seq_meth[i*flac_bin_size + len(seq_meth) - flanking_size: (i+1) * flac_bin_size + len(seq_meth) - flanking_size], threshold)
            if not is_template:
                m_p, m_n = m_n, m_p
            me_df_t[index][i+10] = m_p
            me_df_nt[index][i+10] = m_n
            if m_n != None and m_p != None and i < flac_up_num_binsoverlap:
                flac_up_avg_p[i] = update_mean(flac_up_avg_p[i], flac_up_bins_sum[i], m_p, flac_bin_size)
                flac_up_avg_n[i] = update_mean(flac_up_avg_n[i], flac_up_bins_sum[i], m_n, flac_bin_size)
                flac_up_bins_sum[i] += flac_bin_size
    if exp_prcnt > 0:
        np.save('me_df_t_' + str(exp_prcnt) + 'high'+'.npy', me_df_t)
        np.save('me_df_nt_' + str(exp_prcnt) + 'high'+ '.npy', me_df_nt)
    else:
        np.save('me_df_t_' + str(exp_prcnt) + 'low'+'.npy', me_df_t)
        np.save('me_df_nt_' + str(exp_prcnt) + 'low'+ '.npy', me_df_nt)

    return genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n

#Count_base, flac_num_bin_overlap, all weights are one, method5
def get_gene_meth_count_based3(meth_count_seq, unmeth_count_seq,  genes_df,  bin_num, flanking_size = 2000):
    print('5')
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
        is_template = row['strand'] == '+'
        prev_gene_end = -1
        next_gene_start = sys.maxsize

        if index > 0 and genes_df.iloc[index-1]['chr'] == row['chr']:
            prev_gene_end = genes_df.iloc[index-1]['end']
        if index < len(genes_df)-1 and genes_df.iloc[index+1]['chr'] == row['chr']:
            next_gene_start = genes_df.iloc[index+1]['start']

        seq_meth = meth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
        seq_unmeth = unmeth_count_seq[row['chr']][row['start'] - flanking_size: row['end'] + flanking_size]
        flac_bin_size = int(flanking_size/bin_num)
        gene_bin_size = int((row['end'] - row['start'])/bin_num)
        flac_down_num_binsoverlap = min((row['start'] - prev_gene_end)/flac_bin_size, bin_num)
        flac_up_num_binsoverlap = min((next_gene_start - row['end'])/flac_bin_size, bin_num)


        if not is_template:
            seq_meth = seq_meth[::-1]
            seq_unmeth = seq_unmeth[::-1]
            flac_down_num_binsoverlap, flac_up_num_binsoverlap = flac_up_num_binsoverlap, flac_down_num_binsoverlap

        for i in range(bin_num):
            m_p, m_n = get_meth_percentage_count_based(seq_meth[i*flac_bin_size: (i+1) * flac_bin_size], seq_unmeth[i*flac_bin_size: (i+1) * flac_bin_size])
            if not is_template:
                m_p, m_n = m_n, m_p
            if m_n != None and m_p != None and i > bin_num - flac_down_num_binsoverlap - 1:
                flac_down_avg_p[i] = update_mean(flac_down_avg_p[i], flac_down_bins_sum[i], m_p, 1)
                flac_down_avg_n[i] = update_mean(flac_down_avg_n[i], flac_down_bins_sum[i], m_n, 1)
                flac_down_bins_sum[i] += 1

            m_p, m_n = get_meth_percentage_count_based(seq_meth[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size], seq_unmeth[i*gene_bin_size + flanking_size: (i+1) * gene_bin_size + flanking_size])
            if not is_template:
                m_p, m_n = m_n, m_p
            if m_n != None and m_p != None:
                genes_avg_p[i] = update_mean(genes_avg_p[i], gene_bins_sum[i], m_p, 1)
                genes_avg_n[i] = update_mean(genes_avg_n[i], gene_bins_sum[i], m_n, 1)
                gene_bins_sum[i] += 1

            m_p, m_n = get_meth_percentage_count_based(seq_meth[i*flac_bin_size + len(seq_meth) - flanking_size: (i+1) * flac_bin_size + len(seq_meth) - flanking_size], seq_unmeth[i*flac_bin_size + len(seq_unmeth) - flanking_size: (i+1) * flac_bin_size + len(seq_unmeth) - flanking_size])
            if not is_template:
                m_p, m_n = m_n, m_p
            if m_n != None and m_p != None and i < flac_up_num_binsoverlap:
                flac_up_avg_p[i] = update_mean(flac_up_avg_p[i], flac_up_bins_sum[i], m_p, 1)
                flac_up_avg_n[i] = update_mean(flac_up_avg_n[i], flac_up_bins_sum[i], m_n, 1)
                flac_up_bins_sum[i] += 1

    return genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n

def template_methylation(meth_seq, genes_df, threshold=0.1):
    average_p = 0
    average_n = 0
    sum_weight = 0
    for index, row in genes_df.iterrows():
        meth_stat = meth_seq[row['chr']][row['start'] - 1: row['end'] -1]
        gene_size = row['end'] - row['start']
        is_template = row['strand'] == '+'
        m_p, m_n = get_meth_percentage(meth_stat, threshold)
        if m_p != None and m_n != None:

            if not is_template:
                m_p, m_n = m_n, m_p
            average_p = update_mean(average_p, sum_weight, m_p, gene_size)
            average_n = update_mean(average_n, sum_weight, m_n, gene_size)
            sum_weight += gene_size

    return average_p, average_n

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

def avg_count(sum_count, bin_num):
    avgs = np.zeros(bin_num, dtype=np.double)
    for i in range(bin_num):
        avgs[i] = sum_count[0][i] / (sum_count[0][i] + sum_count[1][i])
    return avgs

def update_mean(mean, sum_weights, new_value, weight):
    return ((mean * sum_weights) + (new_value * weight)) / (sum_weights + weight)

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
