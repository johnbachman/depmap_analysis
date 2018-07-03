import csv
import argparse as ap
import pandas as pd
import depmap_network_functions as dnf
from time import time
import pdb

# Temp fix because indra doesn't import when script is run from terminal =
import sys
sys.path.append('~/repos/indra/')

# 1. no geneset -> use corr from full DepMap data, no filtering needed
# 2. geneset, not strict -> interaction has to contain at least one gene from
#    the loaded list
# 3. geneset loaded, strict -> each correlation can only contain genes from
#    the loaded data set


def read_gene_set_file(gf, data):
    gset = []
    with open(gf, 'rt') as f:
        for g in f.readlines():
            gn = g.upper().strip()
            if gn in data:
                gset.append(gn)
    return gset


def main(args):

    # Open data file
    data = pd.read_csv(args.ceres_file, index_col=0, header=0)
    data = data.T

    # pdb.set_trace()  # Check if 'data' is ok; check args.geneset_file args.strict 
    if args.geneset_file:
        # Read gene set to look at
        gene_set = read_gene_set_file(gf=args.geneset_file, data=data)
        # pdb.set_trace()  # Check if there is a 'gene_set'

    # 1. no loaded gene list OR 2. loaded gene list but not strict -> data.corr
    if not args.geneset_file or (args.geneset_file and not args.strict):
        # Calculate the full correlations, or load from cached
        if args.recalc:
            corr = data.corr()
            corr.to_hdf('correlations.h5', 'correlations')
        else:
            corr = pd.read_hdf('correlations.h5', 'correlations')
        # pdb.set_trace()  # Check corr
        # No gene set file, leave 'corr' intact and unstack
        if not args.geneset_file:
            fcorr_list = corr.unstack()

        # Gene set file present: filter and unstack
        elif args.geneset_file and not args.strict:
            fcorr_list = corr[gene_set].unstack()
            # pdb.set_trace()  # check fcorr_list

    # 3. Strict: both genes in interaction must be from loaded set
    elif args.geneset_file and args.strict:
        fcorr_list = data[gene_set].corr().unstack()

    flarge_corr = fcorr_list[fcorr_list != 1.0]
    flarge_corr = flarge_corr[flarge_corr.abs() > args.ll]
    if args.ul < 1.0:
        flarge_corr = flarge_corr[flarge_corr.abs() < args.ul]
    fsort_corrs = flarge_corr.abs().sort_values(ascending=False)
    # pdb.set_trace()  # Check fsort_corrs

    # Find out if the HGNC pairs are connected and if they are how
    uniq_pairs = set()
    with open(args.outbasename+'_all.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        dir_conn_pairs = []
        unexplained = []
        for pair in fsort_corrs.items():
            (id1, id2), score = pair
            if frozenset([id1, id2, score]) not in uniq_pairs:
                uniq_pairs.add(frozenset([id1, id2, score]))
                wrtr.writerow([id1, id2, score])
                if dnf.are_connected(id1, id2):
                    dir_conn_pairs.append([id1, id2, score,
                                           dnf.connection_type(id1, id2)])
                else:
                    unexplained.append([id1, id2, score])

    with open(args.outbasename+'_connections.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(dir_conn_pairs)

    with open(args.outbasename+'_unexplained.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(unexplained)


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument('-f', '--ceres-file', required=True,
                        help='Ceres correlation score in csv format')
    parser.add_argument('-g', '--geneset-file',
                        help='Filter to interactions with gene set data file.')
    parser.add_argument('-s', '--strict', action='store_true',
                        help='With \'-s\', the correlations are restricted to '
                             'only be between loaded gene set. If no gene set '
                             'is loaded, this option has no effect.')
    parser.add_argument('-stf', '--statement-file',
                        help='With \'-stf\' use file with INDRA statements '
                             'instead of quering from the web clients.')
    parser.add_argument('-o', '--outbasename', default=str(int(time())),
                        help='Base name for outfiles. Default: UTC timestamp')
    parser.add_argument('-r', '--recalc', action='store_true',
                        help='With \'-r\', recalculate full gene-gene '
                             'correlation of data set.')
    parser.add_argument('-ll', type=float, default=0.5,
                        help='Lower limit CERES correlation score filter.')
    parser.add_argument('-ul', type=float, default=1.0,
                        help='Upper limit CERES correlation score filter.')

    a = parser.parse_args()
    main(a)
