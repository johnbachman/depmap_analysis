import sys
import json
import boto3
import logging
import networkx as nx
from networkx import NodeNotFound
from subprocess import call
from datetime import datetime
from flask import Flask, request, abort, Response, redirect, url_for
from indra_db.util import dump_sif
import depmap_network_functions as dnf
from util.io_functions import _pickle_open, _dump_it_to_pickle

logger = logging.getLogger('INDRA GDE API')

# Need way to read/update latest sif dump of pa-statements:
# see indra_db/indra_db/util/dump_sif.py


def _todays_date():
    return datetime.now().strftime('%Y%m%d')


class IndraNetwork:
    """Handle searches and graph output of the INDRA DB network"""
    def __init__(self, indra_graph):
        self.nx_graph_repr = indra_graph
        self.nodes = self.nx_graph_repr.nodes
        self.edges = self.nx_graph_repr.edges

    def find_shortest_path(self, source, target, weight=None):
        try:
            return nx.shortest_path(self.nx_graph_repr, source, target, weight)
        except NodeNotFound:
            return []

    def has_path(self, source, target):
        return nx.has_path(self.nx_graph_repr, source, target)


def dump_indra_db(path='.'):
    base_name = 'db_dump_' + _todays_date()
    if path is not '.':
        path_base_name = path + base_name
    else:
        path_base_name = base_name
    stmts_file = path_base_name + '.pkl'
    dataframe_file = path_base_name + '_dataframe.pkl'
    csv_file = path_base_name + '.csv'
    files = ' '.join((stmts_file, dataframe_file, csv_file))
    cmd = 'python ' + dump_sif.__file__ + ' ' + files
    logger.info('Executing subprocess: %s' % cmd)
    try:
        retcode = call(cmd, shell=True)
        if retcode < 0:
            logger.warning('Script was terminated by signal: ' + str(-retcode))
        else:
            logger.info('Script finished ' + str(retcode))
    except OSError as e:
        logger.error('Script failed: ' + e.strerror)

    return stmts_file, dataframe_file, csv_file

# Need way to create directed graph of sif dump: add function to depmap
# network functions


def load_indra_graph(graph_path, update=False):
    if update:
        stmts_file, dataframe_file, csv_file = dump_indra_db()
        indra_graph = dnf.nx_directed_graph_from_sif_dataframe(dataframe_file)
        logging.info('Dumping latest indra db snapshot to pickle')
        _dump_it_to_pickle(graph_path, indra_graph)
    else:
        indra_graph = _pickle_open(graph_path)
    return indra_graph

# Need way to receive user queries for network search

# Need way to return results
