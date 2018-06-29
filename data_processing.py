from __future__ import print_function
import sys
import itertools
from collections import defaultdict
import numpy as np
import pandas as pd
import sklearn

CONSTRUCTONS = [ 'F9','GP1BB','HBB','HBG1','HNF4A_P2','IRF4','IRF6','LDLR','MSMB',
                 'MYC_rs6983267','PKLR','SORT1','TERT_GBM','TERT_HEK293T','ZFAND3']

def empty_dataframe(data):
    return pd.DataFrame(columns=data.keys())

def select_data_by(data, key, values):
    if values == 'all':
        return data
    else:
        results = []
        for value in values:
            results.append(data[ data[key] == value ])
        return pd.concat(results)

def reject_data_by(data, key, values):
    if values == 'all':
        return empty_dataframe(data)
    else:
        results = []
        for value in np.unique(data[key]):
            if value not in values:
                results.append(data[ data[key] == value ])
        return pd.concat(results)

def select_features(data, features):
    columns = service_column_names(data) + features
    return data[columns]

def service_column_names(data):
    result = ['SNV','construction','value','confidence']
    if 'block_num' in data.keys():
        result.append('block_num')
    return result

def feature_names(data):
    take_always_features = service_column_names(data)
    return list(filter(lambda feature: feature not in take_always_features, data.keys()))

def filter_features(data, condition):
    features_to_add = list(filter(condition, feature_names(data)))
    columns = service_column_names(data) + list(np.unique(features_to_add))
    return data[columns]

def feature_type(feature_name):
    return feature_name.rsplit(':')[-1]

def feature_names_by_type(data, feature_types):
    if feature_types == 'all':
        return feature_names(data)
    else:
        return list(filter(lambda f: feature_type(f) in feature_types, feature_names(data)))

def extract_features_only(data):
    return data[feature_names(data)]

def split_features_and_target(data):
    return (extract_features_only(data), data['confidence'])

def unrelated_constructions(construction):
    if construction.startswith('TERT'):
        return [c for c in CONSTRUCTONS if not c.startswith('TERT')]
    else:
        return [c for c in CONSTRUCTONS if c != construction]

def split_by_construction(data, construction):
    test_data = select_data_by(data, 'construction', [construction])
    train_data = select_data_by(data, 'construction', unrelated_constructions(construction))
    return (train_data, test_data)

def split_by_block(data, block_num):
    test_data = select_data_by(data, 'block_num', [block_num])
    train_data = reject_data_by(data, 'block_num', [block_num])
    return (train_data, test_data)

def filter_features_by_type(data, feature_types):
    return select_features(data, feature_names_by_type(data, feature_types))

