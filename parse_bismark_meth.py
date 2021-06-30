import pandas as pd
import numpy as np

def read_meth(in_addrss, out_addrss):
    methylations = pd.read_csv(in_addrss, sep = '\t')
    methylations.columns = ['chro', 'position', 'strand', 'meth', 'unmeth', 'context', 'tri']
    methylations['meth_lvl'] = methylations.meth.divide(methylations.meth + methylations.unmeth)
    methylations = methylations.replace(np.nan, '-', regex = True)
    methylations.to_csv(out_addrss, sep = '\t')
