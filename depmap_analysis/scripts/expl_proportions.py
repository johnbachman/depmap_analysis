import argparse
import pandas as pd
import matplotlib.pyplot as plt
from os import path
from pathlib import Path

from depmap_analysis.util.io_functions import is_dir_path, pickle_open

# Parameters to care about:
# 1. Graph type
# 2. SD ranges
# 3. Type of explanations


def thousands(n):
    """Turn an int to a string of its value per 1000

    Parameters
    ----------
    n : int

    Returns
    -------
    str
    """
    if n < 1000:
        return str(n)
    else:
        return str(n // 1000) + 'k'


parser = argparse.ArgumentParser()
parser.add_argument('--title', required=True,
                    help='Title for data (will also be user got plot output '
                         'name)')
total_col = 'total checked'
default_cols = ['explained (excl sr)',
                'common parent',
                'explained set',
                'complex or direct',
                'x intermediate']
parser.add_argument('--columns', nargs='+',
                    default=default_cols,
                    help=f'Specify columns to plot. Default: {default_cols}')
parser.add_argument('--explainer-dir', type=is_dir_path(), required=True,
                    help='The explainer files live here. No other pickle '
                         'files should be present.')
parser.add_argument('--labels', nargs='+',
                    help='Legend labels on plot (corresponding to column '
                         'names). Default: column names')
parser.add_argument('--outdir',
                    help='Directory where to put the saved figure. Default '
                         'is same directory as stats file.')
parser.add_argument('--show-plot', action='store_true',
                    help='With this flag active, the generated plots will ' +
                         'be shown as well as saved')


args = parser.parse_args()
expl_dir = Path(args.explainer_dir)
outdir = Path(args.outdir) if args.outdir else expl_dir.joinpath('output')
if not outdir.is_dir():
    outdir.mkdir(parents=True)
data_title = args.title
labels = args.columns
labels = labels if len(labels) > 0 else default_cols
legend_labels = args.labels if args.labels else labels
print('Using legend labels: %s' % ' '.join(legend_labels))

# Store explainers by their graph type
expl_by_type = {'pybel': [],
                'signed': [],
                'unsigned': []}
for explainer_file in expl_dir.glob('*.pkl'):
    expl = pickle_open(explainer_file)
    expl_by_type[expl.script_settings['graph_type']].append(expl)

# Per graph type, extract what the old code has
for graph_type, list_of_explainers in expl_by_type.items():
    if len(list_of_explainers) == 0:
        print(f'Skipping graph type {graph_type}')
        continue
    stats_norm = pd.DataFrame(columns=['range', 'filter_w_count'] + labels)

    for explainer in list_of_explainers:
        sumd = explainer.get_summary()
        N = sumd['total checked']
        data = {k: v/N for k, v in sumd.items() if k in labels}
        lo, hi = explainer.sd_range
        data['range'] = f'{lo}-{hi} SD' if hi else f'{lo}+ SD'
        data['filter_w_count'] = data['range'] + '\n' + thousands(N)
        stats_norm = stats_norm.append(other=pd.DataFrame(data=data,
                                                          index=[0]),
                                       sort=False)
    stats_norm.sort_values('range', inplace=True)

    stats_norm.plot(x='filter_w_count',
                    y=labels,
                    legend=legend_labels,
                    kind='bar',
                    # logy=ylog,
                    title=f'{data_title}, {graph_type.capitalize()}',
                    stacked=False)
    # plt.xticks(rotation=270)
    plt.ylabel('Explained fraction')
    plt.savefig(outdir.joinpath(f'{data_title}_{graph_type}.png'))
    plt.show()

    stats_norm.plot(x='filter_w_count',
                    y=labels,
                    legend=legend_labels,
                    kind='bar',
                    logy=True,
                    title=f'{data_title}, {graph_type.capitalize()} (ylog)',
                    stacked=False)
    # plt.xticks(rotation=270)
    plt.ylabel('Explained fracation')
    plt.savefig(outdir.joinpath(f'{data_title}_{graph_type}_ylog.png'))
    plt.show()
