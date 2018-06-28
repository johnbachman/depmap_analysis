import pandas as pd

data = pd.read_csv('portal-Avana-2018-06-21.csv', index_col=0, header=0)
data = data.transpose()

# The correlations take a long time (> 1 hr) to calculate, so cache them
# and don't recalculate unless desired. They also take some 10 GB in ram 
# to process so be prepared for temporary slowdowns
recalculate = False
if recalculate:
    corr = data.corr()
    corr.to_hdf('correlations.h5', 'correlations')
else:  # Opens previously cached data
    corr = pd.read_hdf('correlations.h5', 'correlations')

#labels_to_drop = get_redundant_pairs(corr)
#au_corr = corr_list.drop(labels=labels_to_drop).sort_values(ascending=False)

corr_list = corr.unstack()  # unstack to df 
large_corr = corr_list[corr_list != 1.0] 
large_corr = large_corr[large_corr.abs() > 0.5]
sort_corrs = large_corr.abs().sort_values(ascending=False)

with open('prior_genes.txt', 'rt') as f:
    prior_genes = [line.strip() for line in f.readlines()]

metab_genes = []
with open('metabolic_genes.txt', 'rt') as f:
    for line in f.readlines():
        gene_name = line.strip().upper()
        if gene_name in data:
            metab_genes.append(gene_name)

prior_corrs = large_corr[metab_genes]  # 
metab_data = data[metab_genes]
metab_corr = metab_data.corr()
mcorr_list = metab_corr.unstack()
mlarge_corr = mcorr_list[mcorr_list != 1.0]
mlarge_corr = mlarge_corr[mlarge_corr.abs() > 0.5]
msort_corrs = mlarge_corr.abs().sort_values(ascending=False)
# metab_data.mean().sort_values(ascending=False)
