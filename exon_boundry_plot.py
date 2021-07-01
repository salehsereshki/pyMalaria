
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyMalaria.constants as constants



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
            if exon[i] == -1 * constants.NON_METH_TAG:
                exon_n[i] = (0 - meth_mean) / meth_std
            else:
                exon_n[i] = (np.abs(exon[i]) - meth_mean) / meth_std
        elif exon[i] > 0:
            exon_p_isc[i] = True
            if exon[i] == constants.NON_METH_TAG:
                exon_p[i] = (0 - meth_mean) / meth_std
            else:
                exon_p[i] = (exon[i] - meth_mean) / meth_std
    return exon_p, exon_p_isc, exon_n, exon_n_isc

def get_Cs_mean_std(exon):
    Cs = []
    for i in range(len(exon)):
        if exon[i] != 0:
            if exon[i] == constants.NON_METH_TAG or exon[i] == -1 * constants.NON_METH_TAG:
                Cs.append(0)
            else:
                Cs.append(abs(exon[i]))

    if len(Cs) == 0:
        raise Exception('No C', 'There is No Cs in the exon interval')
    meth_mean = np.mean(Cs)
    meth_std = np.std(Cs)
    return meth_mean, meth_std



def plot_exon_boundry_5(organism_name, annot_df, meth_seq, boundry=150, from_file = False):
    if not from_file:
        exon5_p, exon5_n, exon3_p, exon3_n = get_exon_boundry_meth_strs(annot_df, meth_seq, boundry_size = boundry)
        save(exon5_p, exon5_n, exon3_p, exon3_n, organism_name)
    else:
        _exons_names = ['exon5_p', 'exon5_n', 'exon3_p', 'exon3_n']
        exon5_p = np.load('./saved_data/'+organism_name+_exons_names[0]+'.npy')
        exon5_n = np.load('./saved_data/'+organism_name+_exons_names[1]+'.npy')
    l_regressions = []

    _m, _b = np.polyfit(range(0, boundry), exon5_p[:boundry], 1)
    l_regressions.append([_m, _b, 'blue', range(0, boundry)])
    _m, _b = np.polyfit(range(0, boundry), exon5_p[boundry:], 1)
    l_regressions.append([_m, _b, 'blue', range(boundry, 2*boundry)])
    _m, _b = np.polyfit(range(0, boundry), exon5_n[:boundry], 1)
    l_regressions.append([_m, _b, 'red', range(0, boundry)])
    _m, _b = np.polyfit(range(0, boundry), exon5_n[boundry:], 1)
    l_regressions.append([_m, _b, 'red', range(boundry, 2*boundry)])


    fig, ax = plt.subplots()
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = True
    plt.grid(False)
    plt.style.use('seaborn')
    plt.box(False)
    ticks = [0, 50, 100, 150, 200, 250, 300]
    lables = ['-150', '-100', '-50', '0', '50', '100', '150']
    plt.xticks(ticks, lables)
    ax.yaxis.tick_left()
    ax.plot(range(0, len(exon5_p)), exon5_p, color='blue', alpha=0.6, linewidth = 0.6)
    ax.plot(range(0, len(exon5_n)), exon5_n, color='red', alpha=0.6, linewidth = 0.6)
    for line in l_regressions:
        ax.plot(line[3], line[0] * line[3]+ line[1], color=line[2], linewidth = 2.0, solid_capstyle='butt')
    ax.axhline(y=0.0, color='black', linestyle='-')
    ax.axvline(x=0.0, color='black', linestyle='-')
    #plt.title('5\' boundary')
    #ax.axvline(x=150.0, color='black', linestyle='-', linewidth=0.2)
    plt.rcParams['axes.facecolor'] = 'white'
    plt.savefig('boundary_5_' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()

def save(exon5_p, exon5_n, exon3_p, exon3_n, organism_name):
    _exons = [exon5_p, exon5_n, exon3_p, exon3_n]
    _exons_names = ['exon5_p', 'exon5_n', 'exon3_p', 'exon3_n']
    for i in range(len(_exons)):
        try:
            np.save('./saved_data/'+organism_name+_exons_names[i]+'.npy', np.array(_exons[i]))
        except TypeError as te:
            print(i)
            print(type(_exons[i]))

def plot_exon_boundry_3(organism_name, annot_df, meth_seq, boundry=150, from_file = False):
    if not from_file:
        exon5_p, exon5_n, exon3_p, exon3_n = get_exon_boundry_meth_strs(annot_df, meth_seq, boundry_size = boundry)
        save(exon5_p, exon5_n, exon3_p, exon3_n, organism_name)
    else:
        _exons_names = ['exon5_p', 'exon5_n', 'exon3_p', 'exon3_n']
        exon3_p = np.load('./saved_data/'+organism_name+_exons_names[2]+'.npy')
        exon3_n = np.load('./saved_data/'+organism_name+_exons_names[3]+'.npy')

    exon3_p = np.nan_to_num(exon3_p, nan=0)
    l_regressions = []

    _m, _b = np.polyfit(range(0, boundry), exon3_p[:boundry], 1)
    l_regressions.append([_m, _b, 'blue', range(-2 * boundry, -boundry)])
    _m, _b = np.polyfit(range(0, boundry), exon3_p[boundry:], 1)
    l_regressions.append([_m, _b, 'blue', range(-boundry, 0)])
    _m, _b = np.polyfit(range(0, boundry), exon3_n[:boundry], 1)
    l_regressions.append([_m, _b, 'red', range(-2 * boundry, -boundry)])
    _m, _b = np.polyfit(range(0, boundry), exon3_n[boundry:], 1)
    l_regressions.append([_m, _b, 'red', range(-boundry, 0)])

    mpl.rcParams.update(mpl.rcParamsDefault)

    fig, ax = plt.subplots()
    ax.grid(False)
    plt.style.use('seaborn')
    plt.box(False)
    ticks = [-300, -250, -200, -150, -100, -50, 0]
    lables = ['-150', '-100', '-50', '0', '50', '100', '150']
    ax.yaxis.tick_right()
    plt.xticks(ticks, lables)
    ax.plot(range(-len(exon3_p), 0), exon3_p, color='blue',  alpha=0.6, linewidth = 0.6)
    ax.plot(range(-len(exon3_p), 0), exon3_n, color='red',  alpha=0.6, linewidth = 0.6)
    for line in l_regressions:
        ax.plot(line[3], line[0] * range(0, boundry)+ line[1], color=line[2], linewidth=2.0, solid_capstyle='butt')

    ax.axhline(y=0.0, color='black', linestyle='-')
    ax.axvline(x=0.0, color='black', linestyle='-')
    #ax.axvline(x=-150.0, color='black', linestyle='-', linewidth=0.2)
    #plt.title('3\' boundary')
    plt.rcParams['axes.facecolor'] = 'white'
    plt.savefig('boundary_3_' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()
