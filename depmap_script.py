import argparse as ap
import pandas as pd
import depmap_network_functions as dnf

# Temp fix because indra doesn't import when script is run from terminal =
import sys
sys.path.append('~/repos/indra/')


def main(args):
    data = pd.read_csv(args.ceres_file, index_col=0, header=0)
    data = data.T

    if args.recalc:
        corr = data.corr()
        corr.to_hdf('correlations.h5', 'correlations')
    else:
        corr = pd.read_hdf('correlations.h5', 'correlations')

    with open('prior_genes.txt', 'rt') as f:
        prior_genes = [line.strip() for line in f.readlines()]
    metab_genes = []
    with open('metabolic_genes.txt', 'rt') as f:
        for line in f.readlines():
            gene_name = line.strip().upper()
            if gene_name in data:
                metab_genes.append(gene_name)

    corr_list = corr.unstack()
    large_corr = corr_list[corr_list != 1.0]
    large_corr = large_corr[large_corr.abs() > args.ll]
    if args.ul < 1.0:
        large_corr = large_corr[large_corr.abs() < args.ul]
    # sort_corrs = large_corr.abs().sort_values(ascending=False)
    # prior_corrs = large_corr[metab_genes]

    metab_data = data[metab_genes]
    metab_corr = metab_data.corr()
    mcorr_list = metab_corr.unstack()
    mlarge_corr = mcorr_list[mcorr_list != 1.0]
    mlarge_corr = mlarge_corr[mlarge_corr.abs() > args.ll]
    if args.ul < 1.0:
        mlarge_corr = mlarge_corr[mlarge_corr.abs() < args.ul]
    msort_corrs = mlarge_corr.abs().sort_values(ascending=False)

    # Find out if the HGCN pairs are connected and if they are how
    dir_conn_pairs = []
    for pair in msort_corrs.items():
        (id1, id2), score = pair
        if dnf.are_connected(id1, id2):
            dir_conn_pairs.append((id1, id2, dnf.connection_type(id1, id2)))

    with open('stats.dat', 'w') as f:
        f.writelines(dir_conn_pairs)


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument('-f', '--ceres-file', required=True,
                        help='Ceres correlation score in csv format')
    parser.add_argument('-r', '--recalc', action='store_true',
                        help='If True (Default), recalculate correlation '
                             'calculation.')
    parser.add_argument('-ll', type=float, default=0.5,
                        help='Lower limit CERES correlation score filter.')
    parser.add_argument('-ul', type=float, default=1.0,
                        help='Upper limit CERES correlation score filter.')

    args = parser.parse_args()
    main(args)
