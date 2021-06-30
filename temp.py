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

def get_exon_boundry_meth_strs(annot_df, meth_seq, boundry_size = 150):
    exons = []
    for chro in meth_seq.keys():
        exons_df = annot_df[(annot_df['type'] == 'exon') & (annot_df['chr'] == chro)].sort_values(by='start')
        genes_df_chr = annot_df[(annot_df['type'] == 'gene') & (annot_df['chr'] == chro)].sort_values(by='start')
        for index, row  in genes_df_chr.iterrows():
            start_index = binarySearch(exons_df['start'].values ,  row['start'])[0]
            gene_exons = []
            while len(exons_df) > start_index and exons_df.iloc[start_index]['start'] >= row['start'] and exons_df.iloc[start_index]['end'] <= row['end']:
                gene_exons.append((exons_df.iloc[start_index]['start'], exons_df.iloc[start_index]['end']))
                start_index += 1
            exons.append([chro, row['strand'],(row['start'], row['end']), gene_exons])

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
                    interval = meth_seq[chro][max(exon_start, exon_end - boundry_size) : min(exon_end + boundry_size, gene_exons[i+1][0])]
                    midpoint = exon_end - max(exon_start, exon_end - boundry_size)
                    exon5.append((interval, midpoint))
                if i > 0:
                    interval = meth_seq[chro][max(exon_start - boundry_size, gene_exons[i-1][1]):min(exon_start + boundry_size , exon_end)]
                    midpoint = exon_start - max(exon_start - boundry_size, gene_exons[i-1][1])
                    exon3.append((interval, midpoint))

        if gene[1] == '-':
            for i in range(len(gene_exons)):
                exon_start = gene_exons[i][0]
                exon_end = gene_exons[i][1]
                chro = gene[0]
                if i > 0:
                    interval = meth_seq[chro][max(gene_exons[i-1][1], exon_start - boundry_size):min(exon_start + boundry_size, exon_end)]
                    midpoint = exon_start - max(gene_exons[i-1][1], exon_start - boundry_size)
                    exon5.append((interval, midpoint))
                if i < len(gene_exons) - 1:
                    interval = meth_seq[chro][max(exon_end-boundry_size, exon_start):min(exon_end + boundry_size, gene_exons[i+1][0])]
                    midpoint = exon_end - max(exon_end-boundry_size, exon_start)
                    exon3.append((interval,midpoint))

    return exon5, exon3


target = []
for exon in exons:
    if exon[1] == '+':
        for j in range(len(exon[3]) - 1):
            exon_start = exon[3][j][0]
            exon_end = exon[3][j][1]
            if sequences[exon[0]][exon[3][j][1] + 2] == 'G' and exon[3][j+1][0] - exon_end > boundry_size and exon_end - exon_start > boundry_size:
                target = exon
                break
    else:
        for j in range(1, len(exon[3])):
            exon_start = exon[3][j][0]
            exon_end = exon[3][j][1]
            if sequences[exon[0]][exon[3][j][0] - 2] == 'G' and exon_start - exon[3][j-1][1] > boundry_size and exon_end - exon_start > boundry_size:
                target = exon
                break

