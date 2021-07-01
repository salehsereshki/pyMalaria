import matplotlib.pyplot as plt
import pandas as pd


def gene_meth_percentage(meth_seq, genes_df,  bin_num, threshold = 0.1):
    genes_meth_p = []
    genes_meth_n = []
    flac_up_p = []
    flac_up_n = []
    flac_down_p = []
    flac_down_n = []
    for index, row in genes_df.iterrows():
        start = row['start']
        end = row['end']
        strand = row['strand']
        fp_start = end
        fp_end = end + 500
        fn_start = start - 500
        fn_end = start
        if strand == '-':
            fp_start = start - 500
            fp_end = start
            fn_start = end
            fn_end = end + 500
        bin_size_f = int(500 / bin_num)
        bin_size = int((end - start) / bin_num)
        meths_p = []
        meths_n = []
        for i in range(bin_num):
            s_i = int(i * bin_size) + start
            e_i = int((i+1) * bin_size) + start
            m_p, m_n = get_meth_percentage(meth_seq[row['chr']][s_i: e_i], threshold)
            if m_p != None and m_n != None:
                meths_p.append(m_p * bin_size)
                meths_n.append(m_n * bin_size)
        meth_up_p = []
        meth_up_n = []
        for i in range(fp_start, fp_end, int(bin_size_f)):
            if i + bin_size_f <= fp_end:
                m_p, m_n = get_meth_percentage(meth_seq[row['chr']][i : i + bin_size_f], threshold)
                if m_p != None and m_n != None:
                    meth_up_p.append(m_p)
                    meth_up_n.append(m_n)
        meth_down_p = []
        meth_down_n = []
        for i in range(fn_start, fn_end, int(bin_size_f)):
            if i + bin_size_f <= fn_end:
                m_p, m_n = get_meth_percentage(meth_seq[row['chr']][i : i + bin_size_f], threshold)
                if m_p != None and m_n != None:
                    meth_down_p.append(m_p)
                    meth_down_n.append(m_n)
        genes_meth_p.append(meths_p)
        genes_meth_n.append(meths_n)
        flac_up_p.append(meth_up_p)
        flac_up_n.append(meth_up_n)
        flac_down_p.append(meth_down_p)
        flac_down_n.append(meth_down_n)
    return genes_meth_p, genes_meth_n, flac_up_p, flac_up_n, flac_down_p, flac_down_n

def get_meth_percentage(meth_stat, threshold):
    countCs_p = 0
    countCs_n = 0
    countMethCs_p = 0
    countMethCs_n = 0
    for i in range(len(meth_stat)):
        if float(meth_stat[i]) > 0:
            countCs_p+=1
        elif float(meth_stat[i]) < 0:
            countCs_n+=1
        if float(meth_stat[i]) > threshold:
            countMethCs_p += 1
        elif float(meth_stat[i]) < -1 * threshold:
            countMethCs_n += 1
    if countCs_p != 0 and countCs_n != 0:
        return float(countMethCs_p) / countCs_p, float(countMethCs_n) / countCs_n
    else:
        return None, None

def get_average(meth_profs_df, genes_df, bin_num, is_flanking):
    meth_df = pd.DataFrame(meth_profs_df, columns = ['0', '1', '2', '3', '4'])
    res = []
    if not is_flanking:
        sum_gene_size = 0
        for index, row in genes_df.iterrows():
            sum_gene_size += int(row['end']) - int(row['start'])
        for i in range(len(meth_df.columns)):
            res.append(meth_df[meth_df.columns[i]].sum()/ (sum_gene_size/bin_num))
    else:
        for i in range(len(meth_df.columns)):
            res.append(meth_df[meth_df.columns[i]].mean())
    return res



def plot_gene_body_meth(organism_name, meth_seq, genes_df, bin_num):
    genes_meth_p, genes_meth_n, flac_up_p, flac_up_n, flac_down_p, flac_down_n = gene_meth_percentage(meth_seq,
                                                                                                              genes_df,
                                                                                                              bin_num)
    # Eech list is consist of all the genes methylation level in each bin. For Gene Body each the vector of size bin_num showing the #mec/#C * bin_size for that gene
    # For Flanking Regions the numbers showing the #meC/#C
    avg_meth_gene_p = get_average(genes_meth_p, genes_df, bin_num, False)
    avg_meth_gene_n = get_average(genes_meth_n, genes_df, bin_num, False)
    avg_meth_up_p = get_average(flac_up_p, genes_df, bin_num, True)
    avg_meth_up_n = get_average(flac_up_n, genes_df, bin_num, True)
    avg_meth_down_p = get_average(flac_down_p, genes_df, bin_num, True)
    avg_meth_down_n = get_average(flac_down_n, genes_df, bin_num, True)

    final_p = avg_meth_down_p + avg_meth_gene_p + avg_meth_up_p
    final_n = avg_meth_down_n + avg_meth_gene_n + avg_meth_up_n

    # yticks = [0]
    # ylabels = ['']
    # plt.yticks(yticks, ylabels)
    plt.tick_params(left=False, labelleft=False)
    plt.box(False)
    plt.ylabel("$meC/C$")
    ticks = [0, 5, 10, 15]
    labels = ['    5\' flanking region', '         gene body', '     3\' flanking region', '']
    plt.xticks(ticks, labels, horizontalalignment='left', verticalalignment='bottom', multialignment='center')
    # plt.tick_params(axis='x', colors='black', direction='out', length=10, width=1,  pad = 4)
    plt.grid(False)
    plt.style.use('seaborn')
    plt.plot(range(0, 15), final_p, color='blue', linewidth=4.0)
    plt.plot(range(0, 15), final_n, color='red', linewidth=4.0)
    plt.axhline(y=0.0, color='black', linestyle='-')
    plt.axvline(x=0.0, color='black', linestyle='-')
    plt.rcParams['axes.facecolor'] = 'white'
    plt.savefig('genebody_' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()
