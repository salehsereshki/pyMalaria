import pyMalaria.methods as methods
import numpy as np
import matplotlib.pyplot as plt


def rgb_to_hex(rgb):
    res = '%02x%02x%02x' % rgb
    return '#'+res

def plot_meth_percent_genome_areas(organism_name, annot_df, sequences, methylations, coverage_threshold, thrshold = 0.1, from_file = False):
    if from_file:
        try:
            res_meth_pct = np.load('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_meth_pct.npy')
            res_pie_pct = np.load('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_pie_pct.npy')
        except (FileNotFoundError, IOError):
            res_meth_pct, res_pie_pct = methods.get_genome_regions_percentage(annot_df, sequences, methylations, coverage_threshold, thrshold = thrshold)
            np.save('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_meth_pct.npy', np.asarray(res_meth_pct))
            np.save('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_pie_pct.npy', np.asarray(res_pie_pct))
    else:
        res_meth_pct, res_pie_pct = methods.get_genome_regions_percentage(annot_df, sequences, methylations, coverage_threshold, thrshold = thrshold)
        np.save('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_meth_pct.npy', np.asarray(res_meth_pct))
        np.save('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_pie_pct.npy', np.asarray(res_pie_pct))

    res_meth_pct = res_meth_pct[0:3]

    fig, ax = plt.subplots(1, 1)
    plt.rcdefaults()
    font = {'size': 14}

    plt.rc('font', **font)
    xlabel = ['exons', 'introns', 'intergenic \n regions']
    #xlabel = ['exons', 'introns', 'intergenic \n regions', 'flanking', 'non flanking']
    x = np.arange(len(res_meth_pct))
    plt.xticks(x + 0.2, xlabel)
    plt.ylabel('% of methylated cytosines')
    plt.tick_params(axis='x', colors='black', direction='out', length=1, width=1,  pad=4, right=True)
    ax.bar(x, res_meth_pct, width=0.4,  align='edge', color=[rgb_to_hex((163, 65, 67)).upper(), rgb_to_hex((67, 94, 156)).upper(), rgb_to_hex((94, 150, 88)).upper()])
    plt.savefig('compartment_meth_proportion' + str(organism_name) + '_thr_'+str(thrshold) +  '.jpg', dpi=2000)
    plt.close()


    plt.rcdefaults()
    font = {'size': 20}

    plt.rc('font', **font)
    labels = ['exons', 'introns', 'intergenic \n regions']
    explode = (0.05, 0.05, 0)
    fig1, ax1 = plt.subplots()
    #colors=[rgb_to_hex((163, 65, 67)).upper(), rgb_to_hex((67, 94, 156)).upper(), rgb_to_hex((94, 150, 88)).upper()]
    #colors = [colors[0], colors[2], colors[1]]
    #res_pie_pct = [res_pie_pct[0], res_pie_pct[2], res_pie_pct[1]]
    #labels = [labels[0], labels[2], labels[1]]
    #ax1.pie(res_pie_pct, explode=explode, labels=labels, autopct='%1.1f%%',
    #    shadow=False, rotatelabels = 270, textprops=dict(color="k"), colors=colors)
    wedges, lbls, txts = ax1.pie(res_pie_pct, explode=explode, labels=labels, autopct='%1.1f%%',labeldistance = 1.3,
        shadow=False, startangle=90, textprops=dict(color="k"), colors=[rgb_to_hex((163, 65, 67)).upper(), rgb_to_hex((67, 94, 156)).upper(), rgb_to_hex((94, 150, 88)).upper()])

    for label in lbls:
        label.set_horizontalalignment('center')
    ax1.axis('equal')
    plt.tight_layout()
    plt.savefig('compartment_meth_repartition' + str(organism_name)+'_thr_'+str(thrshold) + '.jpg', dpi=2000)
    plt.close()
