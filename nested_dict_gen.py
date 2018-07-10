import depmap_network_functions as dnf
from indra.tools.assemble_corpus import load_statements as ac_load_stmts
import argparse as ap
import logging
import pickle
from time import time
from collections import defaultdict
import itertools as itt

logger = logging.info('NestedDictGenerator')

# This script should generate a nested dict strucutre from statements:
# dict(hgnc_id1) -> [dict(hgnc_id2)] | -> [list of connection types]


def main(args):

    # Load statements to create dict structure for
    stmts = ac_load_stmts(args.statements_in)

    stmt_dicts = defaultdict(dict)
    for st in stmts:
        # NOTE1: Agents can be more than two and be only one too.
        # NOTE2: Pair can show up multiple times when connection types differ
        # Hence: Only skip if pair+connection type already exists (frozenset?)

        # Get agent names as list
        agent_names = list(dnf.agent_name_set(st))
        if len(agent_names) > 1:
            # Only connection type per statement
            connection = st.to_json()['type']
            if connection:
                # Permuation: ignore order (i.e. ignore subject/object)
                for agent, other_agent in itt.permutations(agent_names, r=2):
                    try:
                        stmt_dicts[agent][other_agent].add(connection)
                    except KeyError:  # If pair does not exist yet
                        stmt_dicts[agent][other_agent] = {connection}
        else:
            continue
    with open(args.outbasename+'.pkl', 'wb') as pf:
        pickle.dump(obj=stmt_dicts, file=pf)


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument('-o', '--outbasename', default=str(int(time())),
                        help='Base name for outfiles. Default: UTC timestamp')
    parser.add_argument('-sti', '--statements-in', required=True,
                        help='Loads a pickle file to use instead of quering a '
                             'database.')
    parser.add_argument('--verbose', '-v', action='count')
    a = parser.parse_args()
    main(a)
