import pandas as pd
import matplotlib.pyplot as plt
import pdb

file_list = []
next_file = input('path to stats file: ')

while next_file:
    file_list.append(next_file)
    next_file = input('path to stats file: ')

top = input('Cut data at line: ')

data_title = input('Title for data: ')

stats_df = pd.DataFrame()
for n, file in enumerate(file_list):
    current_df = pd.read_csv(file, header=0).head(top) if top else \
        pd.read_csv(file, header=0)
    if n == 0:
        stats_df['filter'] = current_df['filter']
    name = file.split('/')[-1][5:-16]
    stats_df[name] = \
        current_df['direct'] / current_df['total_corr']

stats_df.sort_values('filter', inplace=True)

col_list = stats_df.drop(columns=['filter']).columns.values.tolist()

stats_df.plot(x='filter', y=col_list, legend=col_list, kind='bar',
              title=data_title, stacked=False)

plt.savefig(data_title+'.png')
plt.show()
