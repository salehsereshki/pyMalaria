import pyMalaria.methods as methods
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import os


def plot_exon_boundry_5(organism_name, annot_df, meth_seq, boundry=150, from_file = False):
    if not from_file:
        exon5_p, exon5_n, exon3_p, exon3_n = methods.get_exon_boundry_meth_strs(annot_df, meth_seq, boundry_size = boundry)
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
        exon5_p, exon5_n, exon3_p, exon3_n = methods.get_exon_boundry_meth_strs(annot_df, meth_seq, boundry_size = boundry)
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
