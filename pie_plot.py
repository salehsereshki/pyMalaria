import matplotlib.pyplot as plt

def get_meth_context_percentage(methylations, coverage_thrshld):
    temp = methylations[methylations.meth + methylations.unmeth >= coverage_thrshld]
    temp = temp[temp.meth/ (temp.meth + temp.unmeth) > 0.1]
    temp_cg = temp[temp.context == 'CG']
    temp_chg = temp[temp.context == 'CHG']
    temp_chh = temp[temp.context == 'CHH']
    print('cg', 'chg', 'chh')
    print(len(temp_cg) * 100 / len(temp), len(temp_chg) * 100 / len(temp), len(temp_chh) * 100 / len(temp))
    return([len(temp_cg) * 100 / len(temp), len(temp_chg) * 100 / len(temp), len(temp_chh) * 100 / len(temp)])


def rgb_to_hex(rgb):
    res = '%02x%02x%02x' % rgb
    return '#'+res


def plot_meth_context_percentage(organism_name, methylations, coverage_threshold):

    sizes = get_meth_context_percentage(methylations, coverage_threshold)
    labels = ['CG', 'CHG', 'CHH']
    explode = (0.1, 0.1, 0)
    fig1, ax1 = plt.subplots()
    plt.rcdefaults()
    font = {'size': 14}

    plt.rc('font', **font)
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90, textprops=dict(color="k"), colors=[rgb_to_hex((163, 65, 67)).upper(), rgb_to_hex((67, 94, 156)).upper(), rgb_to_hex((94, 150, 88)).upper()])
    ax1.axis('equal')
    plt.tight_layout()
    plt.savefig('context_meth' + str(organism_name) + '.jpg', dpi=2000)
    plt.close()
