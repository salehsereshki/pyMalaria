
import numpy as np
import matplotlib.pyplot as plt
import input_parser as input_parser


def get_genome_regions_percentage(annot_df, sequences, methylations, coverage_threshold, thrshold = 0.1):

    ## Iterating on the Methylations rows. For the Cytosines which has coverage. The area status of the Cytosine is calculated.
    ##the number of Cytosines with coverage in genomic areas are calculated.
    ##the number of methylated Cytosines with coverage in the genomic areas are calculated.
    ##Divide the number of methylated on each region either by the total number of methylated Cytosines or by the total number of cytosines in that area.
    genes_df = input_parser.subset_genes(annot_df)
    genes_seq = input_parser.make_gene_string(genes_df, sequences)
    exons_df = input_parser.subset_exons(annot_df)
    exons_seq = input_parser.make_exon_string(exons_df, sequences)
    gene_plus_flanking = input_parser.make_gene_plus_flanking_string(genes_df, sequences)

    count_mc_intergenic = 0
    count_c_intergenic = 0
    count_mc_intron = 0
    count_c_intron = 0
    count_mc_exon = 0
    count_c_exon = 0

    count_mc_intergenic_flanking = 0
    count_c_intergenic_flanking = 0

    count_mc_intergenic_fargene = 0
    count_c_intergenic_fargene = 0

    for index, row in methylations.iterrows():
        mc = float(row['meth'])
        unmc = float(row['unmeth'])

        if mc + unmc > coverage_threshold:
            mlvl = mc / (mc + unmc)
            is_gene = genes_seq[row['chr']][row['position'] - 1] == '-1' or genes_seq[row['chr']][row['position'] - 1] == '1'
            is_exon = exons_seq[row['chr']][row['position'] - 1] == '-1' or exons_seq[row['chr']][row['position'] - 1] == '1'
            is_gene_plus_flanking = gene_plus_flanking[row['chr']][row['position'] - 1] == '-1' or gene_plus_flanking[row['chr']][row['position'] - 1] == '1'

            if is_gene and is_exon:
                count_c_exon += 1
                if mlvl > thrshold:
                    count_mc_exon += 1
            elif is_gene and not is_exon:
                count_c_intron += 1
                if mlvl > thrshold:
                    count_mc_intron += 1
            if not is_gene:
                count_c_intergenic += 1
                if mlvl > thrshold:
                    count_mc_intergenic += 1
                if is_gene_plus_flanking:
                    count_c_intergenic_flanking += 1
                    if mlvl > thrshold:
                        count_mc_intergenic_flanking += 1
                else:
                    count_c_intergenic_fargene += 1
                    if mlvl > thrshold:
                        count_mc_intergenic_fargene += 1


    res_meth_pct = [count_mc_exon * 100 / count_c_exon, count_mc_intron * 100 / count_c_intron, count_mc_intergenic * 100 / count_c_intergenic,
                    count_mc_intergenic_flanking * 100 / count_c_intergenic_flanking, count_mc_intergenic_fargene * 100 / count_c_intergenic_fargene
                    ]
    res_pie_pct = [count_mc_exon/ (count_mc_intergenic + count_mc_exon+ count_mc_intron) ,
                   count_mc_intron/ (count_mc_intergenic + count_mc_exon+ count_mc_intron) ,
                   count_mc_intergenic/ (count_mc_intergenic + count_mc_exon+ count_mc_intron) ]
    return res_meth_pct, res_pie_pct


def rgb_to_hex(rgb):
    res = '%02x%02x%02x' % rgb
    return '#'+res

def plot_meth_percent_genome_areas(organism_name, annot_df, sequences, methylations, coverage_threshold, thrshold = 0.1, from_file = False):
    if from_file:
        try:
            res_meth_pct = np.load('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_meth_pct.npy')
            res_pie_pct = np.load('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_pie_pct.npy')
        except (FileNotFoundError, IOError):
            res_meth_pct, res_pie_pct = get_genome_regions_percentage(annot_df, sequences, methylations, coverage_threshold, thrshold = thrshold)
            np.save('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_meth_pct.npy', np.asarray(res_meth_pct))
            np.save('./saved_data/'+organism_name+'_thr_'+str(thrshold)+'_comp_res_pie_pct.npy', np.asarray(res_pie_pct))
    else:
        res_meth_pct, res_pie_pct = get_genome_regions_percentage(annot_df, sequences, methylations, coverage_threshold, thrshold = thrshold)
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
