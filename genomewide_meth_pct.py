import pandas as pd
import sys

def get_genomewide_meth(address):
    methylations = pd.read_table(address)
    methylations.columns = ['chr', 'position', 'strand', 'meth', 'unmeth', 'context', 'three']
    methylation_subset = methylations[methylations['meth'] + methylations['unmeth'] > 3]
    meths = methylation_subset['meth']/ (methylation_subset['meth'] + methylation_subset['unmeth'])
    len(meths)
    print(address, meths.mean() * 100)

#get_genomewide_meth(sys.argv[1])
