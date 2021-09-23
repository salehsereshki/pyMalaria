import input_parser as IP

def get_most_expressed_genes(annot_df, rna_address, percentage):
    exp_df = IP.get_expression(rna_address, annot_df)
    if percentage > 0: #gives the highest expression genes
        exp_df = exp_df[len(exp_df) - int(len(exp_df) * percentage): len(exp_df)]
        exp_df = exp_df.reset_index(drop=True)
        exp_df = exp_df.sort_values(['chr', 'start'], ascending=(True, True))
        exp_df = exp_df.reset_index(drop=True)
        return exp_df
    else:
        percentage = -1*percentage
        exp_df = exp_df[:int(len(exp_df)*percentage)]
        exp_df = exp_df.reset_index(drop=True)
        exp_df = exp_df.sort_values(['chr', 'start'], ascending=(True, True))
        exp_df = exp_df.reset_index(drop=True)
        return exp_df
