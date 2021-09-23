import configs as configs
import gene_body_meth as GBM
import numpy as np
import matplotlib.pyplot as plt
import sys
import gbm_computations as gbmc
import pandas as pd
import seaborn as sns

#This method is for making the genes and flanking region counts dataframe, the columns are
#       ['mp1', 'mp2', 'mp3', 'mp4', 'mp5',
#       'ump1', 'ump2', 'ump3', 'ump4', 'ump5',
#       'mn1', 'mn2', 'mn3', 'mn4', 'mn5',
#       'umn1', 'umn2', 'umn3', 'umn4', 'umn5',
#       'bin_size']
def get_gene_flac_bins_counts_df(meth_count_seq, unmeth_count_seq,  genes_df,  bin_num, flanking_size = 2000):
    genes_counts = np.zeros([len(genes_df), 4 * bin_num + 1], dtype=np.double)
    flac_up_counts = np.zeros([len(genes_df), 4 * bin_num + 1], dtype=np.double)
    flac_down_counts = np.zeros([len(genes_df), 4 * bin_num + 1], dtype=np.double)


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
            flac_down_counts[index][i] = np.sum(sm[sm > 0])
            flac_down_counts[index][bin_num + i] = np.sum(sum[sum > 0])
            flac_down_counts[index][2 * bin_num + i] = abs(np.sum(sm[sm < 0]))
            flac_down_counts[index][3 * bin_num + i] = abs(np.sum(sum[sum < 0]))
            i = i - 1
            bin_start = bin_start - flac_bin_size

        flac_down_counts[index][4 * bin_num] = flac_bin_size

        if row['strand'] == '+':
            bin_start = row['start'] - start_seq_position
        else:
            bin_start = end_seq_position - row['end']
        i = 0
        while i < bin_num:
            sm = seq_meth[bin_start: bin_start + gene_bin_size]
            sum = seq_unmeth[bin_start: bin_start + gene_bin_size]
            genes_counts[index][i] = np.sum(sm[sm > 0])
            genes_counts[index][bin_num + i] = np.sum(sum[sum > 0])
            genes_counts[index][2 * bin_num + i] = abs(np.sum(sm[sm < 0]))
            genes_counts[index][3 * bin_num + i] = abs(np.sum(sum[sum < 0]))
            i = i + 1
            bin_start = bin_start + gene_bin_size
        genes_counts[index][4 * bin_num] = gene_bin_size
        if gene_bin_size == 759:
            print(row, genes_counts[index][2 * bin_num + 1])

        if row['strand'] == '+':
            bin_start = len(seq_meth) - (end_seq_position - row['end'])
        else:
            bin_start = len(seq_meth) - (row['start'] - start_seq_position)
        i = 0
        while i < bin_num and bin_start + flac_bin_size <= len(seq_meth):
            sm = seq_meth[bin_start: bin_start + flac_bin_size]
            sum = seq_unmeth[bin_start: bin_start + flac_bin_size]
            flac_up_counts[index][i] = np.sum(sm[sm > 0])
            flac_up_counts[index][bin_num + i] = np.sum(sum[sum > 0])
            flac_up_counts[index][2 * bin_num + i] = abs(np.sum(sm[sm < 0]))
            flac_up_counts[index][3 * bin_num + i] = abs(np.sum(sum[sum < 0]))
            i = i + 1
            bin_start = bin_start + flac_bin_size
        flac_up_counts[index][4 * bin_num] = flac_bin_size


    columns = ['mp1', 'mp2', 'mp3', 'mp4', 'mp5',
           'ump1', 'ump2', 'ump3', 'ump4', 'ump5',
           'mn1', 'mn2', 'mn3', 'mn4', 'mn5',
           'umn1', 'umn2', 'umn3', 'umn4', 'umn5', 'bin_size']

    return pd.DataFrame(flac_down_counts, columns=columns),\
           pd.DataFrame(genes_counts, columns=columns),\
           pd.DataFrame(flac_up_counts, columns=columns)


#This method finds the correlation between the methylation level of reads in each bin and the size of the bin.
def find_coefficient(bins_counts):
    columns = list(bins_counts.columns)

    for cl in range(5):
        col1 = columns[cl]
        col2 = columns[cl+5]
        col3 = columns[cl+10]
        col4 = columns[cl+15]
        gm1 = bins_counts[col1] / (bins_counts[col1] + bins_counts[col2])
        x, y = [], []
        for i in range(len(gm1)):
            if not gm1.isna()[i]:
                x.append(gm1[i])
                y.append(bins_counts['bin_size'][i])
        print(cl, '1', np.corrcoef(x, y))
        gm2 = bins_counts[col3] / (bins_counts[col3] + bins_counts[col4])
        x, y = [], []
        for i in range(len(gm2)):
            if not gm2.isna()[i]:
                x.append(gm2[i])
                y.append(bins_counts['bin_size'][i])
        print(cl, '2', np.corrcoef(x, y))


#making the dataframe of methylation level and bin size

def get_mlevel_binsize_df(bins_counts):
    columns = list(bins_counts.columns)

    data_p = []
    data_n = []
    for cl in range(5):
        col1 = columns[cl]
        col2 = columns[cl+5]
        col3 = columns[cl+10]
        col4 = columns[cl+15]
        gm1 = bins_counts[col1] / (bins_counts[col1] + bins_counts[col2])
        x, y = [], []
        for i in range(len(gm1)):
            if not gm1.isna()[i]:
                x.append(gm1[i])
                y.append(bins_counts['bin_size'][i])

        data_p.append(pd.DataFrame({'mlevel':x, 'bin_size':y}))
        gm2 = bins_counts[col3] / (bins_counts[col3] + bins_counts[col4])
        x, y = [], []
        for i in range(len(gm2)):
            if not gm2.isna()[i]:
                x.append(gm2[i])
                y.append(bins_counts['bin_size'][i])
        data_n.append(pd.DataFrame({'mlevel':x, 'bin_size':y}))
    return data_p, data_n

def bin_mlevels(mlevel_binsize_df, bin_number, organism_name, region, strand):
    mlevel_binsize_df['bin_tag'] = pd.cut(mlevel_binsize_df['mlevel'], [float(i)/100 for i in range(-1, 99)])
    res = mlevel_binsize_df.groupby('bin_tag').mean()['bin_size'].fillna(0)
    plt.bar(range(len(res)), res)
    plt.savefig('gene_mlevel_binsize_' +str(bin_number)+'_' + region+'_' +strand +'_' +str(organism_name) + '.jpg', dpi=2000)
    plt.close()

#Ploting the box plot of the methylation level of each gene in each bin, ai/(ai+bi) s

def plot_gene_body_box_plot(bins_counts, organism_name, region):
    columns = list(bins_counts.columns)
    data_p = []
    data_n = []
    for cl in range(5):
        col1 = columns[cl]
        col2 = columns[cl+5]
        col3 = columns[cl+10]
        col4 = columns[cl+15]
        gm1 = (bins_counts[col1] / (bins_counts[col1] + bins_counts[col2])).dropna()
        data_p.append(gm1)
        gm2 = (bins_counts[col3] / (bins_counts[col3] + bins_counts[col4])).dropna()
        data_n.append(gm2)
    plt.boxplot(data_p)
    plt.savefig('gene_boxplot_' + region+ '_p_'+str(organism_name) + '.jpg', dpi=2000)
    plt.close()
    plt.boxplot(data_n)
    plt.savefig('gene_boxplot_' + region+ '_n_'+str(organism_name) + '.jpg', dpi=2000)
    plt.close()

    return data_p, data_n


def plot_gene_body_meth(organism_name, meth_seq, genes_df, bin_num, mode, threshold = 0.1, flanking_size = 500, exp_filter=0):
    if mode == 1:
        genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n = gbmc.get_gene_meth_count_based_average_bin(meth_seq[0], meth_seq[1], genes_df,  bin_num, flanking_size=flanking_size)
    if mode == 2:
        genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n = gbmc.get_gene_meth_count_based1(meth_seq[0], meth_seq[1], genes_df,  bin_num, flanking_size=flanking_size)
    if mode == 3:
        genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n = gbmc.get_gene_meth_count_based2(meth_seq[0], meth_seq[1], genes_df,  bin_num, flanking_size=flanking_size)
    if mode == 4:
        genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n = gbmc.get_gene_meth(meth_seq, genes_df,  bin_num, exp_filter, threshold=threshold, flanking_size=flanking_size)
    if mode == 5:
        genes_avg_p, genes_avg_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n = gbmc.get_gene_meth_count_based3(meth_seq[0], meth_seq[1], genes_df,  bin_num, flanking_size=flanking_size)
    final_p = np.concatenate((flac_down_avg_p , genes_avg_p , flac_up_avg_p))
    final_n = np.concatenate((flac_down_avg_n , genes_avg_n , flac_up_avg_n))
    #print(final_p)
    #print(final_n)
    font_size = 16
    yticks = [i*0.002 for i in range(7)]
    ylabels = [str(i*2) for i in range(7)]
    plt.yticks(yticks, ylabels, fontsize=font_size)
    #plt.tick_params(left=False, labelleft=False)
    plt.box(False)
    plt.ylabel("$meC/C(*10^3)$", fontsize=font_size)
    ticks = [0, bin_num, bin_num * 2, bin_num * 3]
    labels = ['   5\' flanking region', '        gene body', '   3\' flanking region', '']
    plt.xticks(ticks, labels, horizontalalignment='left', fontsize=font_size)
    # plt.tick_params(axis='x', colors='black', direction='out', length=10, width=1,  pad = 4)
    plt.grid(False)
    plt.style.use('seaborn')
    plt.plot(range(0, 3 * bin_num), final_p, color='blue', linewidth=5.0)
    plt.plot(range(0, 3 * bin_num), final_n, color='red', linewidth=5.0)
    plt.axhline(y=0.0, color='black', linestyle='-')
    plt.axvline(x=0.0, color='black', linestyle='-')
    plt.axvline(x=5.0, color='black', linestyle='--')
    plt.axvline(x=10.0, color='black', linestyle='--')
    plt.rcParams['axes.facecolor'] = 'white'
    exp_filter_tag = ''
    if exp_filter != 0:
        if exp_filter > 0:
            exp_filter_tag = str(exp_filter)+'high_expr'
        else:
            exp_filter_tag = str( -1 * exp_filter)+'low_expr'

    if mode != 4:
        plt.savefig('genebody_' + str(organism_name) +'_method' + str(mode)+'_fs'+str(flanking_size)+'.jpg', dpi=2000)
    else:
        plt.savefig('genebody_' + str(organism_name) +'_method' + str(mode)+ '_thrsh'+ str(threshold)+'_fs'+str(flanking_size)+exp_filter_tag+'.jpg', dpi=2000)
    plt.close()








