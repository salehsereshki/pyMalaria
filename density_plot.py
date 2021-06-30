import pyMalaria.methods as methods
import pandas as pd
import matplotlib.pyplot as plt
from os import path
import numpy as np


def check_files(organism_name, threshold, chro_num, count_names):
    for name in count_names:
        if 'meth' not in name:
            if not path.exists('./saved_data/' + organism_name + name + '_chro_' + chro_num + '.npy'):
                return False
        else:
            if not path.exists(
                    './saved_data/' + organism_name + name + '_chro_' + chro_num + '_thrsh_' + threshold + '.npy'):
                return False

    return True


def save(organism_name, threshold, chro_num, _counts, _count_names):
    for i in range(len(_counts)):
        if 'meth' not in _count_names[i]:
            np.save('./saved_data/' + organism_name + _count_names[i] + '_chro_' + chro_num + '.npy',
                    np.array(_counts[i]))
        else:
            np.save('./saved_data/' + organism_name + _count_names[
                i] + '_chro_' + chro_num + '_thrsh_' + threshold + '.npy', np.array(_counts[i]))
    return


def load(organism_name, threshold, chro_num, count_names):
    count_Cs_p = np.load('./saved_data/' + organism_name + count_names[0] + '_chro_' + chro_num + '.npy')
    count_Cs_n = np.load('./saved_data/' + organism_name + count_names[1] + '_chro_' + chro_num + '.npy')
    count_methCs_p = np.load(
        './saved_data/' + organism_name + count_names[2] + '_chro_' + chro_num + '_thrsh_' + threshold + '.npy')
    count_methCs_n = np.load(
        './saved_data/' + organism_name + count_names[3] + '_chro_' + chro_num + '_thrsh_' + threshold + '.npy')
    return count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n


def compute_and_save(organism_name, chro, meth_seq, threshold, chro_num, count_names):
    count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n = methods.get_density_vec(chro, meth_seq, thrshld=float(threshold))
    _counts = [count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n]
    save(organism_name, threshold, chro_num, _counts, count_names)
    return count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n


def plot_density_Cs(organism_name, chro, meth_seq, threshold, chro_num, from_file=False):
    count_names = ['count_Cs_p', 'count_Cs_n', 'count_methCs_p', 'count_methCs_n']
    if not from_file or not check_files(organism_name, str(0), str(chro_num), count_names):
        count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n = compute_and_save(organism_name, chro, meth_seq,
                                                                                  str(threshold), str(chro_num),
                                                                                  count_names)
    else:
        count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n = load(organism_name, str(threshold), str(chro_num),
                                                                      count_names)

    x = range(-50, len(count_Cs_p))
    c_p = pd.DataFrame(count_Cs_p)
    c_n = pd.DataFrame(count_Cs_n)
    c_p.columns = ['counts_p']
    c_n.columns = ['counts_n']
    df_c = pd.concat([c_p, c_n], axis=1)
    df_c['counts_n'] = -df_c['counts_n']
    plt.style.use('seaborn')
    plt.rcParams['axes.facecolor'] = 'white'
    df_c.plot(kind='bar', stacked=True, colormap='bwr', grid=False, legend=False)
    plt.xticks(x[::100], rotation='vertical')
    plt.xlabel('position * 1000 (bp)')
    plt.yticks([-500, -250, 0, 250, 500], ['500', '250', '0', '250', '500'])
    plt.ylabel('$C$ count')
    plt.title('chromosome: ' + str(chro_num))
    plt.axhline(y=0.0, color='black', linestyle='-')
    plt.axvline(x=0.0, color='black', linestyle='-')
    plt.tight_layout()
    plt.savefig(organism_name + '_C_chro_' + str(chro_num), dpi=2000)


#        plt.show()

def plot_density_methCs(organism_name, chro, meth_seq, thrshld, chro_num, from_file=False):
    count_names = ['count_Cs_p', 'count_Cs_n', 'count_methCs_p', 'count_methCs_n']

    if not from_file or not check_files(organism_name, str(thrshld), chro_num, count_names):
        count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n = compute_and_save(organism_name, chro, meth_seq,
                                                                                  str(thrshld), str(chro_num),
                                                                                  count_names)
    else:
        count_Cs_p, count_Cs_n, count_methCs_p, count_methCs_n = load(organism_name, str(thrshld), str(chro_num),
                                                                      count_names, count_names)

    x = range(-50, len(count_methCs_p))
    c_p = pd.DataFrame(count_methCs_p)
    c_n = pd.DataFrame(count_methCs_n)
    c_p.columns = ['counts_p']
    c_n.columns = ['counts_n']
    df_c = pd.concat([c_p, c_n], axis=1)
    df_c['counts_n'] = -df_c['counts_n']
    plt.style.use('seaborn')
    plt.rcParams['axes.facecolor'] = 'white'
    df_c.plot(kind='bar', stacked=True, colormap='bwr', grid=False, legend=False)
    plt.xticks(x[::100], rotation='vertical')
    plt.xlabel('position * 1000 (bp)')
    plt.yticks([-40, -20, 0, 20, 40], ['40', '20', '0', '20', '40'])
    plt.ylabel('$me^5C$ count')
    plt.title('chromosome: ' + str(chro_num))
    plt.axhline(y=0.0, color='black', linestyle='-')
    plt.axvline(x=0.0, color='black', linestyle='-')
    plt.tight_layout()
    plt.savefig(organism_name + '_meth_C' + str(thrshld) + '_chro_' + str(chro_num) + '.jpg', dpi=2000)
    plt.close()


def get_partial_meth_C(chro, meth_seq, start, end):
    cnt = 0
    for i in range(start, end):
        if float(meth_seq[chro][i]) < -0.1:
            cnt += 1
    print(cnt)
