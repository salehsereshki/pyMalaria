import matplotlib.pyplot as plt
import pyMalaria.methods as methods


def rgb_to_hex(rgb):
    res = '%02x%02x%02x' % rgb
    return '#'+res


def plot_meth_context_percentage(organism_name, methylations, coverage_threshold):

    sizes = methods.get_meth_context_percentage(methylations, coverage_threshold)
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
