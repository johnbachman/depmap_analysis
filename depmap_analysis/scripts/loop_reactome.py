import logging
import argparse
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from indra.databases.hgnc_client import get_current_hgnc_id, get_uniprot_id
from depmap_analysis.util.io_functions import file_opener, file_path
from depmap_analysis.network_functions.depmap_network_functions import \
    corr_matrix_to_generator, down_sampl_size


logger = logging.getLogger(__name__)


def match_reactome(z_sc, reactome_dict):
    logger.info('Generating generator')
    corr_iterator = corr_matrix_to_generator(z_sc)
    res = {'agA_hgnc': [], 'agA_up': [], 'agB_hgnc': [], 'agB_up': [],
           'z_sc': [], 'has_pathways': [], 'common_pathways': []}
    logger.info('Looping correlations')
    for (a, b), corr in corr_iterator:
        hgnc_id_a = get_current_hgnc_id(a)
        if isinstance(hgnc_id_a, list):
            ix = 0
            while True:
                try:
                    a_up = get_uniprot_id(hgnc_id_a[ix])
                except IndexError:
                    a_up = None
                    break
                if a_up is None:
                    ix += 1
        else:
            a_up = get_uniprot_id(hgnc_id_a)
        if a_up is None:
            continue

        hgnc_id_b = get_current_hgnc_id(b)
        if isinstance(hgnc_id_b, list):
            ix = 0
            while True:
                try:
                    b_up = get_uniprot_id(hgnc_id_b[ix])
                except IndexError:
                    b_up = None
                    break
                if b_up is None:
                    ix += 1
        else:
            b_up = get_uniprot_id(hgnc_id_b)
        if b_up is None:
            continue

        common_reactome = set(reactome_dict.get(a_up, [])) & \
                          set(reactome_dict.get(b_up, []))
        res['agA_hgnc'].append(a)
        res['agA_up'].append(a_up)
        res['agB_hgnc'].append(b)
        res['agB_up'].append(b_up)
        res['z_sc'].append(corr)
        res['common_pathways'].append(common_reactome)
        res['has_pathways'].append(bool(common_reactome))
    logger.info('Returning results')
    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--z-corr', type=file_path('h5'), required=True,
                        help='The path to the stored correlation matrix as '
                             'a pandas DataFrame')
    parser.add_argument('--reactome', type=file_path('pkl'), required=True,
                        help='The reactome pickle')
    args = parser.parse_args()
    z_sc_file = Path(args.z_corr)
    reactome_file = Path(args.reactome)
    sd_ranges = [('rnd', None), (2, 3), (3, 4), (4, 5), (5, None)]

    # Only need first dict
    reactome_mapping = file_opener(reactome_file)[0]

    # Load corr matrix
    z_sc_full = pd.read_hdf(z_sc_file)
    assert isinstance(z_sc_full, pd.DataFrame)

    all_stats = {'range': [], 'checked': [], 'has_pathways': [],
                 'has_pathways_norm': []}
    data_frames = {}

    for ll, ul in sd_ranges:
        # Filter matrix
        if isinstance(ll, (int, float)) and ll and ul:
            logger.info(f'Filtering correlations to {ll} - {ul} SD')
            z_sc_filtered = z_sc_full[((z_sc_full > ll) & (z_sc_full < ul)) |
                                      ((z_sc_full < -ll) & (z_sc_full > -ul))]
        elif isinstance(ll, (int, float)) and ll and not ul:
            logger.info(f'Filtering correlations to {ll}+ SD')
            z_sc_filtered = z_sc_full[(z_sc_full > ll) | (z_sc_full < -ll)]
        elif isinstance(ll, str):
            rnd_sample = 101
            logger.info(f'Doing a random sample of {rnd_sample}')
            z_sc_filtered = z_sc_full.sample(rnd_sample, axis=0)
            z_sc_filtered = z_sc_filtered.filter(list(z_sc_filtered.index),
                                                 axis=1)
        else:
            raise ValueError('Must have both ll and ul defined'
                             ' or set ll to "rnd"')

        # down sample if too large
        target_pairs = 20000
        n_pairs = z_sc_filtered.notna().sum().sum()
        sample_size = down_sampl_size(n_pairs,
                                      len(z_sc_filtered),
                                      target_pairs,
                                      buffer_factor=1.1)
        while n_pairs > 2*target_pairs:
            logger.info(f'Down sampling DataFrame matrix from {n_pairs}')
            z_sc_filtered = z_sc_filtered.sample(sample_size,
                                                 axis=0)
            z_sc_filtered = z_sc_filtered.filter(list(z_sc_filtered.index),
                                                 axis=1)
            n_pairs = z_sc_filtered.notna().sum().sum()
            sample_size = down_sampl_size(n_pairs, len(z_sc_filtered),
                                          target_pairs, buffer_factor=1.1)
        logger.info(f'Correlation matrix sampled to {n_pairs}')

        # Run match
        results = match_reactome(z_sc=z_sc_filtered,
                                 reactome_dict=reactome_mapping)

        # Create df
        res_df = pd.DataFrame(data=results)

        # Do counts
        range_str = f'{ll}-{ul} SD'
        data_frames[range_str] = res_df
        all_stats['range'].append(range_str)
        all_stats['checked'].append(len(res_df))
        all_stats['has_pathways'].append(res_df['has_pathways'].sum())
        all_stats['has_pathways_norm'].append(
            res_df['has_pathways'].sum()/len(res_df)
        )
        logger.info(f'Finished assembling results for {range_str}')

    # Get a stats df
    all_stats_df = pd.DataFrame(data=all_stats)
    all_stats_df.plot(x='range',
                      y=['has_pathways_norm'],
                      legend=['Has pathways'],
                      kind='bar',
                      title='Test title',
                      stacked=False)
    plt.savefig('reactome_matching.png', format='png')
    plt.show()
    plt.close()
    logger.info(f'Saved plot to {Path("reactome_matching.png").absolute()}')
