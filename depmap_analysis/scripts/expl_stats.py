import argparse
import pandas as pd
import matplotlib.pyplot as plt
from os import path

parser = argparse.ArgumentParser('log plotter')
parser.add_argument('--stats-file', help='path to stats file', required=True)
parser.add_argument('--title',
                    help='Title for data (will also be user got plot output '
                         'name)',
                    required=True)
parser.add_argument('--first-n',
                    help='Cut data after n lines. Useful if there are extra '
                         'notes in the csv file after the actual data',
                    type=int)
parser.add_argument('--columns', nargs='+', default=['total_expl_excl_sr',
                                                     'common_parent',
                                                     'mitochondrial',
                                                     'direct',
                                                     'axb_excl_sr'],
                    help='Specify columns to plot. Default: '
                         '[\'total_expl_excl_sr\', \'common_parent\', '
                         '\'mitochondrial\', \'direct\', \'axb_excl_sr\']')
parser.add_argument('--labels', nargs='+',
                    help='Legend labels on plot (corresponding to column '
                         'names). Default: column names')
parser.add_argument('--outdir',
                    help='Directory where to put the saved figure. Default '
                         'is same directory as stats file.')

args = parser.parse_args()
stats_file = args.stats_file
outdir = args.outdir if args.outdir else path.dirname(stats_file)
data_title = args.title
top = args.first_n
labels = args.columns
labels = labels if len(labels) > 0 else ['total_expl_excl_sr',
                                         'common_parent',
                                         'mitochondrial',
                                         'direct',
                                         'axb_excl_sr']
legend_labels = args.labels if args.labels else labels
print('Using legend labels: %s' % ' '.join(legend_labels))

stats_df = pd.read_csv(stats_file, header=0)

stats_df = stats_df.head(top) if top else stats_df
stats_norm = pd.DataFrame()

stats_df.sort_values('range', inplace=True)
stats_norm['range'] = stats_df['range']

# for col in stats_df.drop(columns=['range', 'low_num', 'high_num']).columns:
for col in stats_df.drop(columns=['range']).columns:
    stats_norm[col] = stats_df[col] / stats_df['total_corr']


def thousands(row):
    if row.total_corr < 1000:
        return str(row.total_corr)
    else:
        return str(row.total_corr // 1000) + 'k'


label_count = stats_df.apply(thousands, axis=1)

stats_norm['filter_w_count'] = stats_norm['range'] + '\n' + label_count

stats_norm.plot(x='filter_w_count',
                y=labels,
                legend=legend_labels,
                kind='bar',
                # logy=ylog,
                title=data_title,
                stacked=False)
# plt.xticks(rotation=270)
plt.ylabel('Explained fracation')
plt.savefig(path.join(outdir, data_title+'.png'))
plt.show()

data_title += ' ylog'
stats_norm.plot(x='filter_w_count',
                y=labels,
                legend=legend_labels,
                kind='bar',
                logy=True,
                title=data_title,
                stacked=False)
# plt.xticks(rotation=270)
plt.ylabel('Explained fracation')
plt.savefig(path.join(outdir, data_title+'.png'))
plt.show()
