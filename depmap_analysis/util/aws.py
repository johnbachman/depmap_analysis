import json
import pickle
from operator import itemgetter
from botocore.exceptions import ClientError

from indra.util.aws import get_s3_file_tree, get_s3_client
from indra_depmap_service.util import DUMPS_BUCKET, DUMPS_PREFIX, logger, \
    NET_BUCKET


def get_latest_sif_s3():
    necc_files = ['belief', 'sif', 'src_counts']
    s3 = get_s3_client(unsigned=False)
    tree = get_s3_file_tree(s3, bucket=DUMPS_BUCKET, prefix=DUMPS_PREFIX,
                            with_dt=True)
    # Find all pickles
    keys = [key for key in tree.gets('key') if key[0].endswith('.pkl')]
    # Sort newest first
    keys.sort(key=lambda t: t[1], reverse=True)
    # Get keys of those pickles
    keys_in_latest_dir = [k[0] for k in keys if
                          any(nfl in k[0] for nfl in necc_files)]
    # Map key to resource
    necc_keys = {}
    for n in necc_files:
        for k in keys_in_latest_dir:
            if n in k:
                # Save and continue to next file in necc_files
                necc_keys[n] = k
                break
    df = load_pickle_from_s3(s3, key=necc_keys['sif'],
                             bucket=DUMPS_BUCKET)
    sev = load_pickle_from_s3(s3, key=necc_keys['src_counts'],
                              bucket=DUMPS_BUCKET)
    bd = read_json_from_s3(s3, key=necc_keys['belief'],
                           bucket=DUMPS_BUCKET)
    return df, sev, bd


def load_pickle_from_s3(s3, key, bucket):
    try:
        res = s3.get_object(Key=key, Bucket=bucket)
        pyobj = pickle.loads(res['Body'].read())
    except Exception as err:
        logger.error('Someting went wrong while loading, reading or '
                     'unpickling the object from s3')
        raise err
    return pyobj


def dump_json_to_s3(name, json_obj, public=False, get_url=False):
    """Set public=True for public read access"""
    s3 = get_s3_client(unsigned=False)
    key = 'indra_network_search/' + name
    options = {'Bucket': DUMPS_BUCKET,
               'Key': key}
    if public:
        options['ACL'] = 'public-read'
    s3.put_object(Body=json.dumps(json_obj), **options)
    if get_url:
        return s3.generate_presigned_url(
            'get_object', Params={'Key': key, 'Bucket': DUMPS_BUCKET})


def read_json_from_s3(s3, key, bucket):
    try:
        res = s3.get_object(Key=key, Bucket=bucket)
        json_obj = json.loads(res['Body'].read().decode())
    except Exception as err:
        logger.error('Someting went wrong while loading or reading the json '
                     'object from s3')
        raise err
    return json_obj


def check_existence_and_date_s3(query_hash, indranet_date=None):
    s3 = get_s3_client(unsigned=False)
    key_prefix = 'indra_network_search/%s' % query_hash
    query_json_key = key_prefix + '_query.json'
    result_json_key = key_prefix + '_result.json'
    exits_dict = {}
    if indranet_date:
        # Check 'LastModified' key in results
        # res_query = s3.head_object(Bucket=SIF_BUCKET, Key=query_json_key)
        # res_results = s3.head_object(Bucket=SIF_BUCKET, Key=result_json_key)
        pass
    else:
        try:
            query_json = s3.head_object(Bucket=DUMPS_BUCKET,
                                        Key=query_json_key)
        except ClientError:
            query_json = ''
        if query_json:
            exits_dict['query_json_key'] = query_json_key
        try:
            result_json = s3.head_object(Bucket=DUMPS_BUCKET,
                                         Key=result_json_key)
        except ClientError:
            result_json = ''
        if result_json:
            exits_dict['result_json_key'] = result_json_key
        return exits_dict

    return {}


def dump_pickle_to_s3(name, indranet_graph_object, prefix=''):
    s3 = get_s3_client(unsigned=False)
    key = prefix + name
    s3.put_object(Bucket=NET_BUCKET, Key=key,
                  Body=pickle.dumps(obj=indranet_graph_object))
