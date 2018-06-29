import numpy as np
import pandas as pd
import sklearn
from sklearn import svm
from sklearn.svm import LinearSVC
from sklearn.datasets import make_classification

data = pd.read_csv('motif_features_combined.tsv', sep='\t')
del data['SNV']
del data['construction']
del data['value']
confidence = data['confidence']
answers = 0 + (confidence > 0.1)
del data['confidence']

error_rates = []
for random_state in [13,17,42,48,93,97,99,121,137,1024]:
    features_train, features_test, answers_train, answers_test = sklearn.model_selection.train_test_split(data, answers, test_size=0.05, random_state=random_state)
    clf = LinearSVC(random_state=0)
    clf.fit(features_train, answers_train)
    error_rate = 1 - clf.score(features_test, answers_test)
    error_rates.append(error_rate)
    print(error_rate)
print(np.mean(error_rates), np.std(error_rates))
