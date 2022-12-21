#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import re
import numpy as np

# def select_features(encodings, labels, selected_features, selected_feature_number):
def select_features(encodings,selected_features, selected_feature_number):
    rows = []
    # print("selected_feature",selected_features)
    for i in range(1, len(selected_features)):
        # searchObj = re.search('f\.(\d+)', selected_features[i][0])
        searchObj = re.search('Pos\.(\d+)', selected_features[i][0])

        if searchObj:
            rows.append(int(searchObj.group(1)))
        # print("searchObj", searchObj.group(1))
    feature_number = selected_feature_number
    if len(rows) < selected_feature_number:
        feature_number = len(rows)
    rows = [0] + rows
    # print("row",rows)
    selected_feature_vectors = np.zeros((len(encodings), feature_number + 2)).astype(str)
    # print("ss",type(encodings),np.shape(encodings))
    encodings = np.array(encodings)

    selected_feature_vectors[:, 0] = encodings[:, 0]
    # labels = np.array(['label'] + labels)
    #
    # labels = labels.tolist()
    labels=[0]*(len(encodings)-1)
    # a=len(labels)
    # b=len(labels1)
    labels.insert(0, 'label')
    labels = np.array(labels)
    #
    selected_feature_vectors[:, 1] =labels
    for i in range(1, feature_number + 1):
        selected_feature_vectors[:, i+1] = encodings[:, rows[i]+1]
    return selected_feature_vectors

