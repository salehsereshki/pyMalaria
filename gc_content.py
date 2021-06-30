import numpy as np
import matplotlib.pyplot as plt
import pyMalaria.methods as methods
from scipy import ndimage


def rgb_to_hex(rgb):
    res = '%02x%02x%02x' % rgb
    return '#'+res


def plot_gc_content(gc_percents):
    fig, ax = plt.subplots(1, 1)
    plt.rcdefaults()
    font = {'size': 14}

    plt.rc('font', **font)
    xlabel = ['PC', 'PV', 'PF']
    #xlabel = ['exons', 'introns', 'intergenic \n regions', 'flanking', 'non flanking']
    x = np.arange(len(gc_percents))
    plt.xticks(x + 0.2, xlabel)
    plt.ylabel('CG content percentage')
    plt.tick_params(axis='x', colors='black', direction='out', length=1, width=1,  pad=4, right=True)
    ax.bar(x, gc_percents, width=0.4,  align='edge', color=[rgb_to_hex((163, 65, 67)).upper(), rgb_to_hex((67, 94, 156)).upper(), rgb_to_hex((94, 150, 88)).upper()])
    plt.savefig('cg_content.jpg', dpi=2000)
    plt.close()

def get_cg_percentage(seq):
    if len(seq) == 0:
        return 0.5
    count_c = 0
    count_g = 0
    count_all = 0
    for bp in seq:
        count_all += 1
        if bp.upper() == 'C':
            count_c += 1
        if bp.upper() == 'G':
            count_g += 1

    return (float(count_c + count_g)) / count_all



def gene_cg_percentage(sequences, genes_df, bin_num):
    genes_cg = []
    flac_up_cg = []
    flac_down_cg = []
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
        cgs = []
        for i in range(bin_num):
            s_i = int(i * bin_size) + start
            e_i = int((i+1) * bin_size) + start
            m_p = get_cg_percentage(sequences[row['chr']][s_i: e_i])
            if m_p != None:
                cgs.append(m_p * bin_size)
        cg_up = []
        for i in range(fp_start, fp_end, int(bin_size_f)):
            if i + bin_size_f <= fp_end:
                m_p= get_cg_percentage(sequences[row['chr']][i: i + bin_size_f])
                if m_p != None:
                    cg_up.append(m_p)
        cg_down = []
        for i in range(fn_start, fn_end, int(bin_size_f)):
            if i + bin_size_f <= fn_end:
                m_p = get_cg_percentage(sequences[row['chr']][i: i + bin_size_f])
                if m_p != None:
                    cg_down.append(m_p)
        genes_cg.append(cgs)
        flac_up_cg.append(cg_up)
        flac_down_cg.append(cg_down)
    return genes_cg, flac_up_cg, flac_down_cg


def get_average(genes_cg, flac_up_cg, flac_down_cg, genes_df, bin_num):
    avg_cg_gene = methods.get_average(genes_cg, genes_df, bin_num, False)
    avg_cg_up = methods.get_average(flac_up_cg, genes_df, bin_num, True)
    avg_cg_down = methods.get_average(flac_down_cg, genes_df, bin_num, True)

    final = avg_cg_down + avg_cg_gene + avg_cg_up

    return final

def plot_cg_content_percentage(sequences, genes_df, bin_num, organism_name):
    genes_cg, flac_up_cg, flac_down_cg = gene_cg_percentage(sequences, genes_df, bin_num)
    final = get_average(genes_cg, flac_up_cg, flac_down_cg, genes_df, bin_num)

    #yticks = [0, 20, 40, 60, 80]
    #ylabels = ['0', '20', '40', '60', '80']
    #plt.yticks(yticks, ylabels)
    plt.tick_params(left=False, labelleft=False)
    plt.box(False)
    plt.ylabel("% CG-content")
    ticks = [0, 5, 10, 15]
    labels = ['    5\' flanking region', '         gene body', '     3\' flanking region', '']
    plt.xticks(ticks, labels, horizontalalignment='left', verticalalignment='bottom', multialignment='center')
    # plt.tick_params(axis='x', colors='black', direction='out', length=10, width=1,  pad = 4)
    plt.grid(False)
    plt.style.use('seaborn')
    plt.plot(range(0, 15), final, color='gray', linewidth=6.0)
    plt.axhline(y=0.0, color='black', linestyle='-')
    plt.axvline(x=0.0, color='black', linestyle='-')
    plt.rcParams['axes.facecolor'] = 'white'
    plt.savefig('cg_content_' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()



def plot_all_together(PC_cg, PV_cg, PF_cg):
    #plt.tick_params(left=False, labelleft=False)
    plt.box(False)
    plt.ylabel("% CG-content")
    ticks = [0, 5, 10, 15]
    labels = ['    5\' flanking region', '         gene body', '     3\' flanking region', '']
    plt.xticks(ticks, labels, horizontalalignment='left', verticalalignment='bottom', multialignment='center')
    # plt.tick_params(axis='x', colors='black', direction='out', length=10, width=1,  pad = 4)
    plt.grid(False)
    plt.style.use('seaborn')
    plt.plot(range(0, 15), PC_cg, label='PC', color='teal', linewidth=4.0)
    plt.plot(range(0, 15), PV_cg, label='PV', color='indigo', linewidth=4.0)
    plt.plot(range(0, 15), PF_cg, label='PF', color='firebrick', linewidth=4.0)
    plt.legend()
    plt.axhline(y=0.0, color='black', linestyle='-')
    plt.axvline(x=0.0, color='black', linestyle='-')
    plt.rcParams['axes.facecolor'] = 'white'
    plt.savefig('cg_content_all' + '.jpg', dpi=2000)
    plt.close()
