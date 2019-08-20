import logging
import unittest

from .api import app


logger = logging.getLogger('api test')


class IndraNetworkServiceTest(unittest.TestCase):

    def setUp(self):
        app.testing = True
        self.app = app.test_client()

    def tearDown(self):
        pass

    def query_test(self):
        """Test a query to empty network"""

        query = {
            'source': 'BRCA1',
            'target': 'BRCA2',
            'stmt_filter': [],
            'edge_hash_blacklist': [],
            'node_filter': ['hgnc', 'fplx'],
            'node_blacklist': [],
            'path_length': False,
            'sign': 'no_sign',
            'weighted': False,
            'bsco': 0.0,
            'curated_db_only': False,
            'fplx_expand': False,
            'k_shortest': 1,
            'two_way': False
        }

        resp = self.app.post('query/submit', json=query)
        resp_json = resp.get_json()['result']
        logger.info('Response: %s' % resp_json)
        assert resp.status_code == 200
        assert resp_json.keys() == {'paths_by_node_count', 'common_targets',
                             'common_parents', 'timeout'}
        assert isinstance(resp_json['paths_by_node_count'], dict)
        assert isinstance(resp_json['common_targets'], list)
        assert isinstance(resp_json['common_parents'], dict)
        assert isinstance(resp_json['timeout'], bool)

    def test_query_endpoint(self):
        """Test if endpoint responds

        NOTE: only tests that the query was received by the API and that it
        was properly formatted.
        """

        resp = self.app.post('query/submit', json={'test': 'api'})
        logger.info('api response: %s' % repr(resp))
        assert resp.status_code == 200
        assert resp.get_json()
        assert resp.get_json().get('result', None)
        assert resp.get_json().get('result') == 'api test passed'
