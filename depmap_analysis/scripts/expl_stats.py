import pandas as pd
import matplotlib.pyplot as plt

stats_file = input('path to stats file: ')
data_title = input('Title for data: ')
top = input('Cut data at line: ')
ylog = bool(input('Log scale for y-axis? (press Enter for linear scale) '))
labels = input('Plot which columns? ').split()
labels = labels if len(labels) > 0 else ['total_expl_excl_sr',
                                         'common_parent',
                                         'mitochondrial',
                                         'direct',
                                         'axb_excl_sr']

legend_labels = input('Legend labels (in order of columns, '
                      'split on whitespace): ').split()
legend_labels = legend_labels if len(legend_labels) > 0 else labels
print('Using legend labels: %s' % ' '.join(legend_labels))

stats_df = pd.read_csv(stats_file, header=0)

stats_df = stats_df.head(top) if top else stats_df
stats_norm = pd.DataFrame()

# stats_df.sort_values('filter', inplace=True)
stats_norm['filter'] = stats_df['filter']

# for col in stats_df.drop(columns=['filter', 'low_num', 'high_num']).columns:
for col in stats_df.drop(columns=['filter']).columns:
    stats_norm[col] = stats_df[col] / stats_df['total_corr']


def thousands(row):
    if row.total_corr < 1000:
        return str(row.total_corr)
    else:
        return str(row.total_corr // 1000) + 'k'


label_count = stats_df.apply(thousands, axis=1)

stats_norm['filter_w_count'] = stats_norm['filter'] + '\n' + label_count

stats_norm.plot(x='filter_w_count',
                y=labels,
                legend=legend_labels,
                kind='bar',
                logy=ylog,
                title=data_title,
                stacked=False)
# plt.xticks(rotation=270)
plt.ylabel('Explained fracation')
plt.savefig(data_title+'.png')
plt.show()
