from __future__ import print_function
import sys
import itertools
from collections import defaultdict
import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
import catboost
import matplotlib.pyplot as plt
from data_processing import *

CONSTRUCTONS = [ 'F9','GP1BB','HBB','HBG1','HNF4A_P2','IRF4','IRF6','LDLR','MSMB',
                 'MYC_rs6983267','PKLR','SORT1','TERT_GBM','TERT_HEK293T','ZFAND3']

def get_rocs(class_labels_and_score_pairs, label = 'ROC', filename=None, plot=True):
    COLORS = ['darkorange', 'red', 'green', 'blue', 'cyan', 'magenta']
    results = {}
    lw = 2
    if plot:
        plt.figure()
    for idx,(class_labels, scores, dataset_label) in enumerate(class_labels_and_score_pairs):
        fpr, tpr, _ = sklearn.metrics.roc_curve(class_labels, scores)
        roc_auc = sklearn.metrics.auc(fpr, tpr)
        results[dataset_label] = roc_auc
        if plot:
            plt.plot(fpr, tpr, color=COLORS[idx], lw=lw, 
                     label=dataset_label + ' ROC curve (area = %0.2f)' % roc_auc)
    if plot:
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(label)
        plt.legend(loc="lower right")
        plt.show()
        if filename:
            plt.savefig(filename)  
    return results


def calc_auc(train_data, test_data,
            plot=True, label = 'ROC', filename = None,
            classifierBuilder=lambda:xgb.XGBClassifier(n_estimators=100, max_depth=3, n_jobs=-1), 
            useRegression=False):
    X_train, y_train = split_features_and_target(train_data)
    X_test, y_test = split_features_and_target(test_data)

    if useRegression:
        y_train_to_fit = y_train
    else:
        y_train_to_fit = 0 + (y_train > 0.1)

    y_train_to_assess = 0 + (y_train > 0.1)
    y_test_to_assess = 0 + (y_test > 0.1)

    classifier = classifierBuilder()
    classifier.fit(X_train, y_train_to_fit)
    
    if useRegression:
        y_score_test = classifier.predict(X_test)
        y_score_train = classifier.predict(X_train)
    else:
        y_score_test = classifier.predict_proba(X_test)[:,1]
        y_score_train = classifier.predict_proba(X_train)[:,1]
    
    roc_results = get_rocs(
        [ (y_test_to_assess, y_score_test, 'test'), 
          (y_train_to_assess, y_score_train, 'train')],
        label=label,
        plot = plot,
        filename = filename)
    return {'roc': roc_results, 'classifier': classifier}


# block_data would be splitted into train_addition and test_block
# train_base is concatenated to train_addition and evaluated on test_block
# evaluation resuls over test_blocks are merged to estimate ROC
def calc_auc_by_blocks(train_base, block_data,
                       plot=True, label = 'ROC', filename = None,
                       classifierBuilder=lambda:xgb.XGBClassifier(n_estimators=100, max_depth=3, n_jobs=-1),
                       useRegression=False):
    y_score_test_total = np.array([])
    y_score_train_total = np.array([])
    y_test_to_assess_total = np.array([])
    y_train_to_assess_total = np.array([])
    for blocknum in np.unique(block_data['block_num']):
        train_data_addition, test_block = split_by_block(block_data, blocknum)
        train_data = pd.concat([train_base, train_data_addition])

        X_train, y_train = split_features_and_target(train_data)
        X_test, y_test = split_features_and_target(test_block)
        
        if useRegression:
            y_train_to_fit = y_train
        else:
            y_train_to_fit = 0 + (y_train > 0.1)
        y_train_to_assess = 0 + (y_train > 0.1)
        y_test_to_assess = 0 + (y_test > 0.1)

        classifier = classifierBuilder()
        classifier.fit(X_train, y_train_to_fit)

        if useRegression:
            y_score_test = classifier.predict(X_test)
            y_score_train = classifier.predict(X_train)
        else:
            y_score_test = classifier.predict_proba(X_test)[:,1]
            y_score_train = classifier.predict_proba(X_train)[:,1]

        y_score_test_total = np.hstack([y_score_test_total, y_score_test])
        y_score_train_total = np.hstack([y_score_train_total, y_score_train])
        y_test_to_assess_total = np.hstack([y_test_to_assess_total, y_test_to_assess])
        y_train_to_assess_total = np.hstack([y_train_to_assess_total, y_train_to_assess])
        
        roc_results = get_rocs(
            [ (y_test_to_assess_total, y_score_test_total, 'test'), 
              (y_train_to_assess_total, y_score_train_total, 'train')],
            label=label, plot = plot, filename = filename)
        return {'roc': roc_results, 'classifier': classifier}

def test_feature_subset(full_data, feature_types, constructions_to_test=CONSTRUCTONS,
                        plot=True,
                        classifierBuilder=lambda:xgb.XGBClassifier(n_estimators=100, max_depth=3, n_jobs=-1),
                        useRegression=False):
    auc_results = {}
    data = filter_features_by_type(full_data, feature_types)
    for construction in constructions_to_test:
        train_data, test_data = split_by_construction(data, construction)
        auc_results[construction] = calc_auc(train_data, test_data, plot=plot, label='ROC for ' + construction)['roc']
    return auc_results

def test_feature_subset_by_blocks(full_data, feature_types, constructions_to_test=CONSTRUCTONS,
                        plot=True,
                        classifierBuilder=lambda:xgb.XGBClassifier(n_estimators=100, max_depth=3, n_jobs=-1),
                        useRegression=False):
    auc_results = {}
    data = filter_features_by_type(full_data, feature_types)
    for construction in constructions_to_test:
        train_data, test_data = split_by_construction(data, construction)
        auc_results[construction] = calc_auc_by_blocks(train_data, test_data, plot=plot, label='block-ROC for ' + construction)['roc']
    return auc_results

def report_results(auc_results, rounding=None):
    print('construction', 'test', 'train', sep='\t')
    for k,v in sorted(auc_results.items()):
        if round:
            print(k, round(v['test'], rounding), round(v['train'], rounding), sep='\t')
        else:
            print(k, v['test'], v['train'], sep='\t')


def stable_feature_importances(train_data, num_rounds=20, train_fraction=0.5,
                               classifierBuilder=lambda:xgb.XGBClassifier(n_estimators=100, max_depth=3, n_jobs=-1),
                               useRegression=False):
    features = feature_names(train_data)
    feature_importances = np.zeros( (len(features), ) )
    for i in range(num_rounds):
        classifier = classifierBuilder()
        X, y = split_features_and_target( train_data.sample(frac=train_fraction) )
        if useRegression:
            classifier.fit(X, 0 + (y > 0.1))
        else:
            classifier.fit(X, y)
        feature_importances += classifier.feature_importances_
        print('.', end='')
    return feature_importances

def sorted_features(feature_names, feature_importances):
    return list(sorted(zip(feature_names,feature_importances), key=lambda kv: -kv[1]))

# class ClassifierSelectingFeatures:
#     def __init__(self, base_classifier_generator, num_features):
#         self.base_classifier_generator = base_classifier_generator
#         self.base_classifier = self.base_classifier_generator()

#     def fit(self, X, y):
#         self.base_classifier.fit(X, y)

#     def predict(self, X):
#         self.base_classifier.predict(X)
    
#     def predict_proba(self, X):
#         self.base_classifier.predict_proba(X)
