import numpy as np
import pandas as pd
import sklearn
from sklearn import svm
from sklearn.svm import LinearSVC
from sklearn.datasets import make_classification
from sklearn import linear_model, decomposition, datasets
from __future__ import print_function
from collections import defaultdict
import matplotlib.pyplot as plt
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_curve, auc
import itertools
import catboost

colors = ['darkorange', 'red', 'green', 'blue', 'cyan', 'magenta']

constructions = [ 'F9','GP1BB','HBB','HBG1','HNF4A_P2','IRF4','IRF6','LDLR','MSMB',
                 'MYC_rs6983267','PKLR','SORT1','TERT_GBM','TERT_HEK293T','ZFAND3']

def remove_unnecessary_columns(data_full):
    data = data_full.copy()
    del data['SNV']
    del data['value']
    del data['construction']
    confidence = data['confidence']
    answers = 0 + (confidence > 0.1)
    del data['confidence']
    return (data, answers)
    
data_full = pd.read_csv('motif_features_combined_fc.tsv', sep='\t')



feature_ranks = {}
classifiers = {}
for construction in constructions:
    X_train, y_train = remove_unnecessary_columns(data_full[data_full['construction'] != construction])
    clf = catboost.CatBoostClassifier().load_model('models/catboost/' + construction + '.cbm')
    feature_importances = clf.get_feature_importance(X_train, y_train)
#     print('===========================\n'+construction)
#     for feature, importance in itertools.islice(sorted(zip(data.keys(), feature_importances), key=lambda pair: -pair[1]), 0, 10):
#         print(feature, round(importance, 3), sep="\t")
    for rank, (feature, importance) in enumerate(sorted(zip(X_train.keys(), feature_importances), key=lambda pair: -pair[1])):
        if feature not in feature_ranks:
            feature_ranks[feature] = []
        feature_ranks[feature].append(rank)

features_to_add = []
for feature, ranks in itertools.islice(sorted(feature_ranks.items(), key= lambda pair: np.mean(pair[1])), 50):
    # print(feature, np.median(ranks), np.mean(ranks), np.std(ranks))
    features_to_add.append(feature)


data_full = data_full.filter(items = ['SNV', 'construction', 'value', 'confidence'] + features_to_add)



auc_results = {}
for construction in constructions:
    X_train, y_train = remove_unnecessary_columns(data_full[data_full['construction'] != construction])
    X_test, y_test = remove_unnecessary_columns(data_full[data_full['construction'] == construction])
    
#     feature_selector_1 = sklearn.feature_selection.VarianceThreshold()
#     feature_selector_1.fit(X_train, y_train)
#     X_train = feature_selector_1.transform(X_train)
#     X_test = feature_selector_1.transform(X_test)
    
#     feature_selector_2 = sklearn.feature_selection.SelectKBest(sklearn.feature_selection.f_classif, k=10)
#     feature_selector_2.fit(X_train, y_train)
#     X_train = feature_selector_2.transform(X_train)
#     X_test = feature_selector_2.transform(X_test)

#     lsvc = LinearSVC(C=1, penalty="l1", dual=False, verbose=1)
#     lsvc.fit(X_train, y_train)    
#     model_select = sklearn.feature_selection.SelectFromModel(lsvc, prefit=True)
#     X_train = model_select.transform(X_train)
#     X_test = model_select.transform(X_test)

    #X_train = pca.transform(X_train)
    #X_test = pca.transform(X_test)    
    
    #X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(data, answers, train_size=0.5,
    #                                                    random_state=rand_state)
    
    classifier = svm.LinearSVC(class_weight='balanced', random_state=13,C=1)
    classifier.fit(X_train, y_train)
    y_score_test = classifier.decision_function(X_test)
    y_score_train = classifier.decision_function(X_train)
    
#     classifier = catboost.CatBoostClassifier(random_seed=17, logging_level='Verbose')
#     classifier.fit(X_train, y_train)
#     #classifier.save_model('models/catboost/'+construction+'.cbm')
#     y_score_test = classifier.predict_proba(X_test)[:,1]
#     y_score_train = classifier.predict_proba(X_train)[:,1]
    
    fpr_test, tpr_test, _ = roc_curve(y_test, y_score_test)
    roc_auc_test = auc(fpr_test, tpr_test)
      
    fpr_train, tpr_train, _ = roc_curve(y_train, y_score_train)
    roc_auc_train = auc(fpr_train, tpr_train)

    lw = 2
    plt.figure()
    plt.plot(fpr_test, tpr_test, color=colors[0],
             lw=lw, label='test ROC curve (area = %0.2f)' % roc_auc_test)
    plt.plot(fpr_train, tpr_train, color=colors[1],
             lw=lw, label='train ROC curve (area = %0.2f)' % roc_auc_train)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC for ' + construction)
    plt.legend(loc="lower right")
    #plt.savefig('roc_curves/' + construction + '.png')
    plt.show()
    #print(construction + ": ", "train -", round(roc_auc_train,3), '; test -', round(roc_auc_test,3))
    auc_results[construction] = {'train': roc_auc_train, 'test': roc_auc_test}

print(auc_results)
