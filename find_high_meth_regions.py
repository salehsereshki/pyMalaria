import exon_boundry_plot as EBP
import pandas as pd
import density_plot as DP

def get_high_meth_intervals(count_meth_p, counts_meth_n, methC_in_kpb = 30):
    intervals = []
    if len(count_meth_p) != len(counts_meth_n):
        print('WARNING: TWO COUNT METHYLATION VECTORS ARE NOT SAME SIZE')
    for interval in range(len(count_meth_p)):
        if counts_meth_n[interval] + count_meth_p[interval] > methC_in_kpb:
            intervals.append((interval*1000, (interval + 1) * 1000, counts_meth_n[interval] + count_meth_p[interval]))
    return intervals


def find_closet_genes(intervals, annot_df, chro_name, meth_seq, threshold = 0.1):
    genes = annot_df[annot_df['chr'] == chro_name]
    genes = genes[genes['type'] == 'gene']
    genes = genes.sort_values('start')
    genes['mid'] = (genes['end'] - genes['start']) / 2 + genes['start']

    df = pd.DataFrame([], columns=['chro', 'geneID', 'strand', 'start', 'end', 'gene meth', 'close interval', 'interval meth in 1 kbp'])

    for i in range(len(intervals)):
        interval_mid = int((intervals[i][1] - intervals[i][0]) / 2 +intervals[i][0])
        index = EBP.binarySearch(list(genes['mid']), interval_mid)[0]
        try:
            cur_gene = genes.iloc[index]
        except IndexError:
            print(intervals[i], len(intervals), index, len(genes))
            continue
        cur_gene_indx = index
        while cur_gene['end'] > intervals[i][0] - 1000 and cur_gene['end'] < intervals[i][1]:
            df = df.append({
                'chro': cur_gene['chr'],
                'geneID': cur_gene['attributes'],
                'strand': cur_gene['strand'],
                'start': cur_gene['start'],
                'end': cur_gene['end'],
                'gene meth': get_gene_meth(cur_gene['chr'], cur_gene['start'], cur_gene['end'], meth_seq, threshold),
                'close interval': (intervals[i][0], intervals[i][1]),
                'interval meth in 1 kbp': intervals[i][2],
            }, ignore_index=True)
            cur_gene_indx -= 1
            if cur_gene_indx < 0:
                break
            cur_gene = genes.iloc[cur_gene_indx]
        cur_gene = genes.iloc[index + 1]
        cur_gene_indx = index + 1
        while cur_gene['start'] < intervals[i][1] + 1000 and cur_gene['start'] < intervals[i][0]:
            df = df.append({
                'chro': cur_gene['chr'],
                'geneID': cur_gene['attributes'],
                'strand': cur_gene['strand'],
                'start': cur_gene['start'],
                'end': cur_gene['end'],
                'gene meth': get_gene_meth(cur_gene['chr'], cur_gene['start'], cur_gene['end'], meth_seq, threshold),
                'close interval': (intervals[i][0], intervals[i][1]),
                'interval meth in 1 kbp': intervals[i][2],
            }, ignore_index=True)
            cur_gene_indx += 1
            if cur_gene_indx == len(genes):
                break
            cur_gene = genes.iloc[cur_gene_indx]

    return df

def get_gene_meth(chr, start, end, meth_seq, threshold = 0.1):
    meth_str = meth_seq[chr][start:end]
    count_meth = 0
    count_C = 0
    for m in meth_str:
        m = float(m)
        if m != 0.0:
            count_C += 1
            if m < -1 * threshold or m > threshold:
                count_meth += 1
    return "{:.2f}".format(float(count_meth) / count_C)


def run(meth_seq, annot_df, organism_name, chromosomes, meth_threshold = 0.1, meth_C_in_kbp = 15):
    df = pd.DataFrame([], columns=['chro', 'geneID', 'strand', 'start', 'end', 'gene meth', 'close interval', 'interval meth in 1 kbp'])
    for indx, chr in enumerate(chromosomes):
        count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n = DP.load(organism_name, meth_threshold, indx+1, DP.count_names)
        intervals = get_high_meth_intervals(count_methCs_p, count_methCs_n, methC_in_kpb = meth_C_in_kbp)
        df = pd.concat([df, find_closet_genes(intervals, annot_df, str(chr), meth_seq, threshold = meth_threshold)])
    df.to_csv('close_genes_'+ organism_name+ '.csv')


