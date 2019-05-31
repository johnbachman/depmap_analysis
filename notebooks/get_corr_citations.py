import csv
from time import sleep
from indra.literature.pubmed_client import get_id_count

def save(citation_counts, filename):
    with open(filename, 'wt') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerows(citation_counts)

if __name__ == '__main__':
    corrs = []
    with open('test_dir_path_all_correlations.csv', 'rt') as f:
        csvreader = csv.reader(f, delimiter=',')
        for row in csvreader:
            corrs.append(row)

    citation_counts = []

    for g1, g2, corr in corrs[18000:]:
        search_query = '%s AND %s' % (g1, g2)
        cite_count = get_id_count(search_query)
        citation_counts.append((g1, g2, corr, str(cite_count)))
        print("%d. %s: %d" % (len(citation_counts), search_query, cite_count))
        if len(citation_counts) % 1000 == 0:
            save(citation_counts, 'citation_counts_%d.csv' % len(citation_counts))
        sleep(0.5)
    save(citation_counts, 'citation_counts_%d.csv' % len(citation_counts))


