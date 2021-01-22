"""Dumps new indra graphs from the latest sif dump"""
import argparse
import depmap_analysis.network_functions.net_functions as nf
from depmap_analysis.util.aws import get_latest_sif_s3, dump_pickle_to_s3, \
    NEW_NETS_PREFIX

__all__ = ['INDRA_MDG', 'INDRA_DG', 'INDRA_SNG', 'INDRA_SEG', 'INDRA_PBSNG',
           'INDRA_PBSEG']

INDRA_MDG = 'indranet_multi_digraph.pkl'
INDRA_DG = 'indranet_dir_graph.pkl'
INDRA_SNG = 'indranet_sign_node_graph.pkl'
INDRA_SEG = 'indranet_sign_edge_graph.pkl'
INDRA_PBSNG = 'indranet_sign_node_pybel.pkl'
INDRA_PBSEG = 'indranet_sign_edge_pybel.pkl'


def dump_new_nets(mdg: bool = False, dg: bool = False, sg: bool = False,
                  spbg: bool = False, verbosity: int = 0,
                  add_mesh_ids: bool = False):
    """Main script function for dumping new networks from latest db dumps"""
    options = dict()

    if add_mesh_ids:
        df, mid = get_latest_sif_s3(get_mesh_ids=True)
        mid_dict = dict()
        for pair in mid:
            mid_dict.setdefault(pair[0], []).append(pair[1])
        options['mesh_id_dict'] = mid_dict
    else:
        df = get_latest_sif_s3()

    options.update({'df': df, 'include_entity_hierarchies': True,
                    'verbosity': verbosity})

    if mdg:
        network = nf.sif_dump_df_to_digraph(graph_type='multi', **options)
        dump_pickle_to_s3(INDRA_MDG, network, prefix=NEW_NETS_PREFIX)
    if dg:
        network = nf.sif_dump_df_to_digraph(**options)
        dump_pickle_to_s3(INDRA_DG, network, prefix=NEW_NETS_PREFIX)
    if sg:
        network, isng = nf.sif_dump_df_to_digraph(graph_type='signed',
                                                  **options)
        dump_pickle_to_s3(INDRA_SEG, network, prefix=NEW_NETS_PREFIX)
        dump_pickle_to_s3(INDRA_SNG, isng, prefix=NEW_NETS_PREFIX)
    if spbg:
        pb_seg, pb_sng = nf.db_dump_to_pybel_sg()
        dump_pickle_to_s3(INDRA_SNG, pb_sng, prefix=NEW_NETS_PREFIX)
        dump_pickle_to_s3(INDRA_SEG, pb_seg, prefix=NEW_NETS_PREFIX)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Dump new networks')
    parser.add_argument('--mdg', help='Dump a new MultiDiGraph',
                        action='store_true', default=False)
    parser.add_argument('--dg', help='Dump a new DiGraph',
                        action='store_true', default=False)
    parser.add_argument('--sg', help='Dump new signed edge and node graphs',
                        action='store_true', default=False)
    parser.add_argument('--pb', help='Dump new PyBel signed edge and node '
                                     'graphs',
                        action='store_true', default=False)
    parser.add_argument('--s3', help='Also upload the new graphs to s3',
                        action='store_true', default=False)
    parser.add_argument('--overwrite', help='Overwrite local files',
                        action='store_true', default=False)
    args = parser.parse_args()
    dump_new_nets(mdg=args.mdg, dg=args.dg, sg=args.sg, spbg=args.pb)
