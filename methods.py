from Bio import SeqIO
import pandas as pd
import numpy as np

NON_METH_TAG = 0.00000001
thrshld = 0.1
def readfasta(address):
    cowpeaRecords = SeqIO.parse(address, "fasta")
    sequences = {}
    for chro in cowpeaRecords:
        sequences[chro.id] = chro.seq
    for i in sequences.keys():
        sequences[i] = sequences[i].upper()
    return sequences
def read_annot(address):
    annot_df = pd.read_table(address, sep='\t', comment='#')
    annot_df.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    return annot_df

def read_methylations(address):
    methylations = pd.read_table(address)
    methylations.columns = ['chr', 'position', 'strand', 'meth', 'unmeth', 'context', 'three']
    return methylations

def check_annot_chro_names(annot_df, sequences):
    annot_chrs = annot_df['chr'].unique()
    count = 0
    for chr in annot_chrs:
        if chr in sequences.keys():
            count += 1
    if len(annot_df['chr'].unique()) == count:
        return True
    else:
        return False

def subset_genes(annot_df):
    genes_df = annot_df[annot_df['type'] == 'gene']
    return genes_df[['chr', 'strand', 'start', 'end']]

def subset_exons(annot_df):
    gene_exon_df = annot_df[(annot_df['type'] == 'exon')]
    return gene_exon_df[['chr', 'strand', 'start', 'end']]

def make_gene_string(genes_df, sequences):
    genes_seq = {}
    for chr in sequences.keys():
        genes_seq[chr] = list(''.zfill(len(sequences[chr])))
    for index, row in genes_df.iterrows():
        gene_marker = '1'
        if row['strand'] == '-':
            gene_marker = '-1'
        for i in range(int(row['start'] - 1), int(row['end'] - 1)):
            genes_seq[row['chr']][i] = gene_marker
    return genes_seq

def make_gene_plus_flanking_string(genes_df, sequences):
    make_gene_plus_flanking = {}
    for chr in sequences.keys():
        make_gene_plus_flanking[chr] = list(''.zfill(len(sequences[chr])))
    for index, row in genes_df.iterrows():
        gene_marker = '1'
        if row['strand'] == '-':
            gene_marker = '-1'
        for i in range(int(row['start']) - 500, int(row['end']) + 500):
            if i < len(make_gene_plus_flanking[row['chr']]) and i > 0:
                make_gene_plus_flanking[row['chr']][i] = gene_marker
    return make_gene_plus_flanking

def make_exon_string(exons_df, sequences):
    exons_seq = {}
    for chr in sequences.keys():
        exons_seq[chr] = list(''.zfill(len(sequences[chr])))
    for index, row in exons_df.iterrows():
        gene_marker = '1'
        if row['strand'] == '-':
            gene_marker = '-1'
        for i in range(int(row['start']), int(row['end'])):
            exons_seq[row['chr']][i] = gene_marker
    return exons_seq

def get_compartment_density(annot_df, sequences, chro):
    sliding_window = 1000
    genes_df = subset_genes(annot_df)
    genes_seq = make_gene_string(genes_df, sequences)

    is_gens = []
    for i in range(0, len(genes_seq[chro]), sliding_window):
        if genes_seq[chro][i: i + sliding_window].count('0') > sliding_window/2 :
            is_gens.append(0)
        else:
            is_gens.append(1)
    return is_gens

def get_genome_dist_str(annot_df):
    genes_df = subset_genes(annot_df)
    exons_df = subset_exons(annot_df)



def get_genome_regions_percentage(annot_df, sequences, methylations, coverage_threshold, thrshold = 0.1):

    ## Iterating on the Methylations rows. For the Cytosines which has coverage. The area status of the Cytosine is calculated.
    ##the number of Cytosines with coverage in genomic areas are calculated.
    ##the number of methylated Cytosines with coverage in the genomic areas are calculated.
    ##Divide the number of methylated on each region either by the total number of methylated Cytosines or by the total number of cytosines in that area.
    genes_df = subset_genes(annot_df)
    genes_seq = make_gene_string(genes_df, sequences)
    exons_df = subset_exons(annot_df)
    exons_seq = make_exon_string(exons_df, sequences)
    gene_plus_flanking = make_gene_plus_flanking_string(genes_df, sequences)

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


def get_index_context(context):
    if context == 'CG':
        return 0
    elif context == 'CHH':
        return 1
    elif context == 'CHG':
        return 2


def get_index(is_gene,  gene_strand, is_cytosine, cytosine_strand, is_methylated):
    digits = [0, 0, 0, 0, 0]
    if is_gene: digits[0] = 1
    if gene_strand == '+': digits[1] = 1
    if is_cytosine: digits[2] = 1
    if cytosine_strand == '+': digits[3] = 1
    if is_methylated: digits[4] = 1
    indx = 16 * digits[0] + 8 * digits[1] + 4 * digits[2] + 2 * digits[3] + digits[4]
    return indx

def get_meth_percentage_gene_non_gene(counts):
    count_gene_meth_cp_gp = counts[get_index(True, '+', True, '+', True)]
    count_gene_unmeth_cp_gp = counts[get_index(True, '+', True, '+', False)]
    count_gene_meth_cn_gn = counts[get_index(True, '-', True, '-', True)]
    count_gene_unmeth_cn_gn = counts[get_index(True, '-', True, '-', False)]
    similar_strand_gene_meth_percent = (count_gene_meth_cp_gp + count_gene_meth_cn_gn) * 100 / (count_gene_meth_cp_gp + count_gene_unmeth_cp_gp + count_gene_meth_cn_gn + count_gene_unmeth_cn_gn)

    count_gene_meth_cp_gn = counts[get_index(True, '-', True, '+', True)]
    count_gene_unmeth_cp_gn = counts[get_index(True, '-', True, '+', False)]
    count_gene_meth_cn_gp = counts[get_index(True, '+', True, '-', True)]
    count_gene_unmeth_cn_gp = counts[get_index(True, '+', True, '-', False)]
    oppose_strand_gene_meth_percent = (count_gene_meth_cp_gn + count_gene_meth_cn_gp) * 100 / (count_gene_meth_cp_gn + count_gene_unmeth_cp_gn + count_gene_meth_cn_gp + count_gene_unmeth_cn_gp)

    non_gene_meth = counts[get_index(False, '+', True, '+', True)] + counts[get_index(False, '+', True, '-', True)]
    non_gene_unmeth = counts[get_index(False, '+', True, '+', False)] + counts[get_index(False, '+', True, '-', False)]

    non_gene_meth_percentage = non_gene_meth * 100 / (non_gene_meth + non_gene_unmeth)
    return similar_strand_gene_meth_percent, oppose_strand_gene_meth_percent, non_gene_meth_percentage

def get_meth_percentage_all_gene_non_gene(cg_counts, chh_counts, chg_counts):
    count_gene_meth_cp_gp = cg_counts[get_index(True, '+', True, '+', True)] + chh_counts[get_index(True, '+', True, '+', True)] + chg_counts[get_index(True, '+', True, '+', True)]
    count_gene_unmeth_cp_gp = cg_counts[get_index(True, '+', True, '+', False)] + chh_counts[get_index(True, '+', True, '+', False)] + chg_counts[get_index(True, '+', True, '+', False)]
    count_gene_meth_cn_gn = cg_counts[get_index(True, '-', True, '-', True)] + chh_counts[get_index(True, '-', True, '-', True)] + chg_counts[get_index(True, '-', True, '-', True)]
    count_gene_unmeth_cn_gn = cg_counts[get_index(True, '-', True, '-', False)] + chh_counts[get_index(True, '-', True, '-', False)] + chg_counts[get_index(True, '-', True, '-', False)]
    similar_strand_gene_meth_percent = (count_gene_meth_cp_gp + count_gene_meth_cn_gn) * 100 / (count_gene_meth_cp_gp + count_gene_unmeth_cp_gp + count_gene_meth_cn_gn + count_gene_unmeth_cn_gn)

    count_gene_meth_cp_gn = cg_counts[get_index(True, '-', True, '+', True)] + chh_counts[get_index(True, '-', True, '+', True)] + chg_counts[get_index(True, '-', True, '+', True)]
    count_gene_unmeth_cp_gn = cg_counts[get_index(True, '-', True, '+', False)] + chh_counts[get_index(True, '-', True, '+', False)] + chg_counts[get_index(True, '-', True, '+', False)]
    count_gene_meth_cn_gp = cg_counts[get_index(True, '+', True, '-', True)] + chh_counts[get_index(True, '+', True, '-', True)] + chg_counts[get_index(True, '+', True, '-', True)]
    count_gene_unmeth_cn_gp = cg_counts[get_index(True, '+', True, '-', False)] + chh_counts[get_index(True, '+', True, '-', False)] + chg_counts[get_index(True, '+', True, '-', False)]
    oppose_strand_gene_meth_percent = (count_gene_meth_cp_gn + count_gene_meth_cn_gp) * 100 / (count_gene_meth_cp_gn + count_gene_unmeth_cp_gn + count_gene_meth_cn_gp + count_gene_unmeth_cn_gp)

    non_gene_meth = cg_counts[get_index(False, '+', True, '+', True)] + cg_counts[get_index(False, '+', True, '-', True)] + chh_counts[get_index(False, '+', True, '+', True)] + chh_counts[get_index(False, '+', True, '-', True)]+ chg_counts[get_index(False, '+', True, '+', True)] + chg_counts[get_index(False, '+', True, '-', True)]
    non_gene_unmeth = cg_counts[get_index(False, '+', True, '+', False)] + cg_counts[get_index(False, '+', True, '-', False)] + chh_counts[get_index(False, '+', True, '+', False)] + chh_counts[get_index(False, '+', True, '-', False)] + chg_counts[get_index(False, '+', True, '+', False)] + chg_counts[get_index(False, '+', True, '-', False)]
    non_gene_meth_percentage = non_gene_meth * 100 / (non_gene_meth + non_gene_unmeth)
    return similar_strand_gene_meth_percent, oppose_strand_gene_meth_percent, non_gene_meth_percentage



def temp_gene_meth_count(gene_seq, meth_seq, context_seq, threshold = 0.1):
    chh_counts = [0] * 32
    cg_counts = [0] * 32
    chg_counts = [0] * 32
    count = 0
    seq_len = 0
    for i in meth_seq.keys():
        seq_len += len(meth_seq[i])
    for chro in meth_seq.keys():
        for i in range(len(meth_seq[chro])):
            count += 1
            if count % 1000000 == 0:
                print(count/1000000, seq_len/1000000)

            is_gene = False
            is_cytosine = False
            gene_strand = '+'
            cytosine_strand = '+'
            is_methylated = False
            if gene_seq[chro][i] != '0':
                is_gene = True
                if gene_seq[chro][i] == '-1':
                    gene_strand = '-'
            if meth_seq[chro][i] != '0':
                is_cytosine = True
                if float(meth_seq[chro][i]) < 0:
                    cytosine_strand = '-'
                if float(meth_seq[chro][i]) > threshold or float(meth_seq[chro][i]) < -1 * threshold:
                    is_methylated = True
            indx = get_index(is_gene,  gene_strand, is_cytosine, cytosine_strand, is_methylated)
            if context_seq[chro][i] == '2':
                chg_counts[indx] = chg_counts[indx] + 1
            elif context_seq[chro][i] == '3':
                chh_counts[indx] = chh_counts[indx] + 1
            else:
                cg_counts[indx] = cg_counts[indx] + 1
    cg_pct_gene_similar, cg_pct_gene_oppose, cg_pct_nongene = get_meth_percentage_gene_non_gene(cg_counts)
    chg_pct_gene_similar, chg_pct_gene_oppose, chg_pct_nongene = get_meth_percentage_gene_non_gene(chg_counts)
    chh_pct_gene_similar, chh_pct_gene_oppose, chh_pct_nongene = get_meth_percentage_gene_non_gene(chh_counts)
    all_pct_gene_similar, all_pct_gene_oppose, all_pct_nongene = get_meth_percentage_all_gene_non_gene(cg_counts, chg_counts, chh_counts)
    print('gene_similar_strand    gene_oppose_strand      non_gene')
    print('cg', cg_pct_gene_similar, cg_pct_gene_oppose, cg_pct_nongene)
    print('chg', chg_pct_gene_similar, chg_pct_gene_oppose, chg_pct_nongene)
    print('chh', chh_pct_gene_similar, chh_pct_gene_oppose, chh_pct_nongene)
    print('all', all_pct_gene_similar, all_pct_gene_oppose, all_pct_nongene)


def get_gene_meth_count(seq_address, meth_address, annot_address, threshold = 0.1):
    sequences = readfasta(seq_address)
    annot_df = read_annot(annot_address)
    methylations = read_methylations(meth_address)
    genes_df = subset_genes(annot_df)
    genes_seq = make_gene_string(genes_df, sequences)

    meth_gene_c = [0, 0, 0]
    unmeth_gene_c = [0, 0, 0]
    meth_nongene_c = [0, 0, 0]
    unmeth_nongene_c = [0, 0, 0]


    thrshld = threshold
    for index, row in methylations.iterrows():
        if genes_seq[row['chr']][int(row['position'])] == '0':
            if int(row['meth']) + int(row['unmeth']) > 0 and float(row['meth']) / (int(row['meth']) + int(row['unmeth'])) > thrshld:
                meth_nongene_c[get_index_context(row['context'])] = meth_nongene_c[get_index_context(row['context'])] + 1
            else:
                unmeth_nongene_c[get_index_context(row['context'])] = unmeth_nongene_c[get_index_context(row['context'])] + 1
        else:
            if int(row['meth']) + int(row['unmeth']) > 0 and float(row['meth']) / (int(row['meth']) + int(row['unmeth'])) > thrshld:
                meth_gene_c[get_index_context(row['context'])] = meth_gene_c[get_index_context(row['context'])] + 1
            else:
                unmeth_gene_c[get_index_context(row['context'])] = unmeth_gene_c[get_index_context(row['context'])] + 1

    all_counts = [meth_gene_c, unmeth_gene_c, meth_nongene_c, unmeth_nongene_c]
    for row in all_counts:
        row.append(row[0]+row[1]+row[2])
    print('\t CG \t CHH \t CHG \t ALL' + '\n')
    print('m_gene' + str(all_counts[0]))
    print('unm_gene' + str(all_counts[1]))
    print('m_ngene' + str(all_counts[2]))
    print('unm_ngene' + str(all_counts[3]))
    gene_meth_percents = []
    nongene_meth_percents = []
    for i in range(len(meth_gene_c)):
        gene_meth_percents.append(float(all_counts[0][i])/(all_counts[0][i] + all_counts[1][i]))
        nongene_meth_percents.append(float(all_counts[2][i])/(all_counts[2][i] + all_counts[3][i]))
    print('gene meth percents : '+ str(gene_meth_percents))
    print('non gene meth percents : ' + str(nongene_meth_percents))
    return(all_counts, gene_meth_percents, nongene_meth_percents)

def make_meth_string(methylations, sequences, coverage_thrshld):
    meth_seq = {}
    #context_seq = {}
    for chr in sequences.keys():
        meth_seq[chr] = list(''.zfill(len(sequences[chr])))
        #context_seq[chr] = list(''.zfill(len(sequences[chr])))
    for index, row in methylations.iterrows():
        meth_c = float(row['meth'])
        unmeth_c = float(row['unmeth'])
        coverage = meth_c + unmeth_c
        strand = row['strand']
        if meth_c > 0 and coverage >= coverage_thrshld:
            if strand == '+':
                meth_seq[row['chr']][row['position'] - 1] = meth_c / (meth_c + unmeth_c)
            else:
                meth_seq[row['chr']][row['position'] - 1] = -meth_c / (meth_c + unmeth_c)
        elif coverage >= coverage_thrshld:
            if strand == '+':
                meth_seq[row['chr']][row['position'] - 1] = NON_METH_TAG
            else:
                meth_seq[row['chr']][row['position'] - 1] = -1 * NON_METH_TAG

        #cntxt = row['context']
        #    cntx_marker = '0'
        #if cntxt == 'CG':
        #    cntx_marker = '1'
        #elif cntxt == 'CHG':
        #    cntx_marker = '2'
        #elif cntxt == 'CHH':
        #    cntx_marker = '3'
        #context_seq[row['chr']][row['position'] - 1] = cntx_marker
    return meth_seq, None

def find_longest_chro(sequences):
    max_lenth = 0
    chro_id = ''
    for chro in sequences.keys():
        if len(sequences[chro]) > max_lenth:
            max_lenth = len(sequences[chro])
            chro_id = chro
    return chro_id

def get_density_vec(chro_name, meth_seq, thrshld):
    sliding_window = 1000
    chro = chro_name
    count_Cs_p = []
    count_Cs_n = []
    count_methCs_p = []
    count_methCs_n = []
    for i in range(0, len(meth_seq[chro]), sliding_window):
        seq = meth_seq[chro]
        meth_c_p = 0
        meth_c_n = 0
        count_c_p = 0
        count_c_n = 0
        for j in range(i, i + sliding_window):
            if j < len(seq):
                if float(seq[j]) > thrshld:
                    meth_c_p += 1
                elif float(seq[j]) < -1 * thrshld:
                    meth_c_n += 1

                if float(seq[j]) > 0:
                    count_c_p += 1
                elif float(seq[j]) < 0:
                    count_c_n += 1

        count_methCs_p.append(meth_c_p)
        count_methCs_n.append(meth_c_n)
        count_Cs_p.append(count_c_p)
        count_Cs_n.append(count_c_n)
    return np.array(count_Cs_p), np.array(count_Cs_n), np.array(count_methCs_p), np.array(count_methCs_n)




def get_meth_context_percentage(methylations, coverage_thrshld):
    temp = methylations[methylations.meth + methylations.unmeth >= coverage_thrshld]
    temp = temp[temp.meth/ (temp.meth + temp.unmeth) > 0.1]
    temp_cg = temp[temp.context == 'CG']
    temp_chg = temp[temp.context == 'CHG']
    temp_chh = temp[temp.context == 'CHH']
    print('cg', 'chg', 'chh')
    print(len(temp_cg) * 100 / len(temp), len(temp_chg) * 100 / len(temp), len(temp_chh) * 100 / len(temp))
    return([len(temp_cg) * 100 / len(temp), len(temp_chg) * 100 / len(temp), len(temp_chh) * 100 / len(temp)])



def get_coverage(methylations):
    meth_sum = methylations['meth'].sum()
    unmeth_sum = methylations['unmeth'].sum()
    print((meth_sum + unmeth_sum)/len(methylations))

#methylations = read_methylations(meth_address2)
#get_coverage(methylations)




def get_meth_percentage(meth_stat):
    countCs_p = 0
    countCs_n = 0
    countMethCs_p = 0
    countMethCs_n = 0
    for i in range(len(meth_stat)):
        if float(meth_stat[i]) > 0:
            countCs_p+=1
        elif float(meth_stat[i]) < 0:
            countCs_n+=1
        if float(meth_stat[i]) > thrshld:
            countMethCs_p += 1
        elif float(meth_stat[i]) < -1 * thrshld:
            countMethCs_n += 1
    if countCs_p != 0 and countCs_n != 0:
        return float(countMethCs_p) / countCs_p, float(countMethCs_n) / countCs_n
    else:
        return None, None

def gene_meth_percentage(meth_seq, genes_df, bin_num):
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
            m_p, m_n = get_meth_percentage(meth_seq[row['chr']][s_i: e_i])
            if m_p != None and m_n != None:
                meths_p.append(m_p * bin_size)
                meths_n.append(m_n * bin_size)
        meth_up_p = []
        meth_up_n = []
        for i in range(fp_start, fp_end, int(bin_size_f)):
            if i + bin_size_f <= fp_end:
                m_p, m_n = get_meth_percentage(meth_seq[row['chr']][i : i + bin_size_f])
                if m_p != None and m_n != None:
                    meth_up_p.append(m_p)
                    meth_up_n.append(m_n)
        meth_down_p = []
        meth_down_n = []
        for i in range(fn_start, fn_end, int(bin_size_f)):
            if i + bin_size_f <= fn_end:
                m_p, m_n = get_meth_percentage(meth_seq[row['chr']][i : i + bin_size_f])
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



def make_exon_boundry_vec(annot_df):
    gene_exon_df = subset_genes_exons(annot_df)
    for index, row in gene_exon_df.itterrows():
        if row['type'] == 'gene':
            gene_exons = []
            ex_i = index + 1
            while gene_exon_df.iloc[ex_i]['type'] == 'exon':
                gene_exons.append(gene_exon_df.iloc[ex_i])
                ex_i += 1
            for i in range(len(gene_exons) - 1):
                start = gene_exons[i]['start']
                end = gene_exons[i]['end']
                flac_up_start = end
                flac_up_end = end + 150
                flac_down_start = start - 150
                flac_down_end = start
                #for j in range()


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
    if mth == NON_METH_TAG:
        return 0
    else:
        return mth

def get_exon_output(exon_meth_start_p, exon_meth_start_n, exon_meth_end_p, exon_meth_end_n, flac_up_avg_p, flac_up_avg_n, flac_down_avg_p, flac_down_avg_n):
    down_p = flac_down_avg_p + exon_meth_start_p
    down_n = flac_down_avg_n + exon_meth_start_n

    up_p = exon_meth_end_p + flac_up_avg_p
    up_n = exon_meth_end_n + flac_up_avg_n

    return down_p, down_n, up_p, up_n


def binarySearch(data, val):
    highIndex = len(data)-1
    lowIndex = 0
    while highIndex > lowIndex:
            index = int((highIndex + lowIndex) / 2)
            sub = data[index]
            if data[lowIndex] == val:
                    return [lowIndex, lowIndex]
            elif sub == val:
                    return [index, index]
            elif data[highIndex] == val:
                    return [highIndex, highIndex]
            elif sub > val:
                    if highIndex == index:
                            return sorted([highIndex, lowIndex])
                    highIndex = index
            else:
                    if lowIndex == index:
                            return sorted([highIndex, lowIndex])
                    lowIndex = index
    return sorted([highIndex, lowIndex])

def get_gene_exons(annot_df, meth_seq):
    exons = []
    for chro in meth_seq.keys():
        exons_df = annot_df[(annot_df['type'] == 'exon') & (annot_df['chr'] == chro)].sort_values(by='start')
        genes_df_chr = annot_df[(annot_df['type'] == 'gene') & (annot_df['chr'] == chro)].sort_values(by='start')
        for index, row in genes_df_chr.iterrows():
            start_index = binarySearch(exons_df['start'].values ,  row['start'])[0]
            gene_exons = []
            while len(exons_df) > start_index and exons_df.iloc[start_index]['start'] >= row['start'] and exons_df.iloc[start_index]['end'] <= row['end']:
                gene_exons.append((exons_df.iloc[start_index]['start'], exons_df.iloc[start_index]['end']))
                start_index += 1
            exons.append([chro, row['strand'], (row['start'], row['end']), gene_exons])

    return exons

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

def get_exon_boundry_meth_strs(annot_df, meth_seq, boundry_size = 150):
    exons = get_gene_exons(annot_df, meth_seq)
    #np.save('exons.npy', np.asarray(exons))
    exon5 = []
    exon3 = []
    for gene in exons:
        gene_exons = gene[3]
        if gene[1] == '+':
            for i in range(len(gene_exons)):
                exon_start = gene_exons[i][0]
                exon_end = gene_exons[i][1]
                chro = gene[0]
                if i < len(gene_exons) - 1:
                    interval = meth_seq[chro][max(exon_start, exon_end - boundry_size): min(exon_end + boundry_size, gene_exons[i+1][0])] # plus one is to make the interval
                    midpoint = exon_end - max(exon_start, exon_end - boundry_size)
                    if (gene_exons[i+1][0] - exon_end > boundry_size and True) and (exon_end - exon_start > boundry_size and True):
                        exon5.append((interval, midpoint))
                if i > 0:
                    interval = meth_seq[chro][max(exon_start - boundry_size, gene_exons[i-1][1]):min(exon_start + boundry_size, exon_end)]
                    midpoint = exon_start - max(exon_start - boundry_size, gene_exons[i-1][1])
                    if (exon_start - gene_exons[i-1][1] > boundry_size and True) and (exon_end - exon_start > boundry_size and True):
                        exon3.append((interval, midpoint))


        if gene[1] == '-':
            for i in range(len(gene_exons)):
                exon_start = gene_exons[i][0]
                exon_end = gene_exons[i][1]
                chro = gene[0]
                if i > 0:
                    interval = meth_seq[chro][max(gene_exons[i-1][1], exon_start - boundry_size):min(exon_start + boundry_size, exon_end)]
                    midpoint = min(exon_start + boundry_size, exon_end) - exon_start - 1
                    if (exon_start - gene_exons[i-1][1] > boundry_size and True) and (exon_end - exon_start > boundry_size and True):
                        exon5.append((interval[::-1], midpoint)) # We need to reverse the exon boundary for the negative strand exons. After this the first of these exon boundary will match to the first of the exon bounaries of the positive strand
                if i < len(gene_exons) - 1:
                    interval = meth_seq[chro][max(exon_end-boundry_size, exon_start):min(exon_end + boundry_size, gene_exons[i+1][0])]
                    midpoint = min(exon_end + boundry_size, gene_exons[i+1][0]) - exon_end - 1
                    if (gene_exons[i+1][0] - exon_end > boundry_size and True) and (exon_end - exon_start > boundry_size and True):
                        exon3.append((interval[::-1], midpoint))


    avg_exon5_p, avg_exon5_n = get_exon_intron_bdr_avg(exon5, boundry_size)
    avg_exon3_p, avg_exon3_n = get_exon_intron_bdr_avg(exon3, boundry_size)
    return avg_exon5_p, avg_exon5_n, avg_exon3_p, avg_exon3_n

def get_Cs_mean_std(exon):
    Cs = []
    for i in range(len(exon)):
        if exon[i] != 0:
            if exon[i] == NON_METH_TAG or exon[i] == -1 * NON_METH_TAG:
                Cs.append(0)
            else:
                Cs.append(abs(exon[i]))

    if len(Cs) == 0:
        raise Exception('No C', 'There is No Cs in the exon interval')
    meth_mean = np.mean(Cs)
    meth_std = np.std(Cs)
    return meth_mean, meth_std

def get_normalized_exon(exon):
    exon = np.asarray(exon, dtype='float')

    exon_p = np.zeros(len(exon))
    exon_n = np.zeros(len(exon))
    exon_p_isc = [False] * len(exon)
    exon_n_isc = [False] * len(exon)

    meth_mean, meth_std = get_Cs_mean_std(exon)

    if meth_std == 0:
        raise Exception('Std zero', 'The std of Cs is zero')

    for i in range(len(exon)):
        if exon[i] < 0:
            exon_n_isc[i] = True
            if exon[i] == -1 * NON_METH_TAG:
                exon_n[i] = (0 - meth_mean) / meth_std
            else:
                exon_n[i] = (np.abs(exon[i]) - meth_mean) / meth_std
        elif exon[i] > 0:
            exon_p_isc[i] = True
            if exon[i] == NON_METH_TAG:
                exon_p[i] = (0 - meth_mean) / meth_std
            else:
                exon_p[i] = (exon[i] - meth_mean) / meth_std
    return exon_p, exon_p_isc, exon_n, exon_n_isc


def get_exon_intron_bdr_avg(meth_strs, boundry_size):
    avgs_p = np.zeros(boundry_size * 2)
    avgs_n = np.zeros(boundry_size * 2)
    sum_p = np.zeros(boundry_size * 2)
    sum_n = np.zeros(boundry_size * 2)
    counts_p = np.zeros(boundry_size * 2)
    counts_n = np.zeros(boundry_size * 2)

    for ex_i in range(len(meth_strs)):
        exon = meth_strs[ex_i][0]
        midpoint = meth_strs[ex_i][1]
        try:
            exon_p, exon_p_isc, exon_n, exon_n_isc = get_normalized_exon(exon)
            for exon_index in range(len(exon)):
                avg_index = exon_index + boundry_size - midpoint
                if exon_p_isc[exon_index]:
                    sum_p[avg_index] += exon_p[exon_index]
                    counts_p[avg_index] += 1
                if exon_n_isc[exon_index]:
                    sum_n[avg_index] += exon_n[exon_index]
                    counts_n[avg_index] += 1

        except Exception as e:
            continue


    for i in range(len(sum_p)):
        if counts_p[i] == 0:
            print("run time warning is because p: ", i)
        avgs_p[i] = sum_p[i] / counts_p[i]
    for i in range(len(sum_n)):
        if counts_n[i] == 0:
            print("run time warning is because n: ", i)
        avgs_n[i] = sum_n[i] / counts_n[i]

    return avgs_p, avgs_n

def compute_gc_content(genes_df, sequences):
    count_c = 0
    count_g = 0
    count_all = 0
    for index, row in genes_df.iterrows():
        for i in range(row['start'] - 1, row['end'] - 1):
            bp = sequences[row['chr']][i]
            count_all += 1
            if bp.upper() == 'C':
                count_c += 1
            if bp.upper() == 'G':
                count_g += 1

    return (float(count_c + count_g)) / count_all

