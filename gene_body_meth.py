import pyMalaria.methods as methods
import matplotlib.pyplot as plt


def plot_gene_body_meth(organism_name, meth_seq, genes_df, bin_num):
    genes_meth_p, genes_meth_n, flac_up_p, flac_up_n, flac_down_p, flac_down_n = methods.gene_meth_percentage(meth_seq,
                                                                                                              genes_df,
                                                                                                              bin_num)

    # Eech list is consist of all the genes methylation level in each bin. For Gene Body each the vector of size bin_num showing the #mec/#C * bin_size for that gene
    # For Flanking Regions the numbers showing the #meC/#C
    avg_meth_gene_p = methods.get_average(genes_meth_p, genes_df, bin_num, False)
    avg_meth_gene_n = methods.get_average(genes_meth_n, genes_df, bin_num, False)
    avg_meth_up_p = methods.get_average(flac_up_p, genes_df, bin_num, True)
    avg_meth_up_n = methods.get_average(flac_up_n, genes_df, bin_num, True)
    avg_meth_down_p = methods.get_average(flac_down_p, genes_df, bin_num, True)
    avg_meth_down_n = methods.get_average(flac_down_n, genes_df, bin_num, True)

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

# ticks = [0, 75, 150, 225, 300]
# labels = ['', '', '0', '75', '150']
# plt.xticks(ticks, labels)
# plt.plot(range(0,300), down_p)
# plt.plot(range(0,300), down_n)
# plt.show()
