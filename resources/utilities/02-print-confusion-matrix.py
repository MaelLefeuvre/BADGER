#!/usr/bin/env python

import sys
from sklearn.metrics import confusion_matrix
import numpy as np
import pandas as pd
def grups_to_degree(pred):
    translator = {
            "Cousins"      : "Third",
            "Third_degree" : "Third",
            "GPGC"         : "Second",
            "Half-siblings": "Second",
            "Second_degree": "Second",
            "Siblings"     : "First",
            "First_degree" : "First",
            "Self"         : "Self",
            "Twins_or_self": "Self",
            "Unrelated"    : "Unrelated"
    }
    return translator[pred]

def make_confusion_matrix(truth_values, predictions):
    return confusion_matrix(truth_values, predictions, labels=np.unique(truth_values))

def print_confusion_matrix(confusion):
    print('   ' + '  '.join([str(n) for n in range(confusion.shape[1])]))
    for rownum in range(confusion.shape[0]):
        print(str(rownum) + '  ' + '  '.join([str(n) for n in confusion[rownum]]))

def pd_print_confusion_matrix(truth_values, predictions):
    labels = np.unique(truth_values)
    matrix = make_confusion_matrix(truth_values, predictions)
    df     = pd.DataFrame(matrix, index=labels, columns=labels)
    #for label in labels:
    #    if label not in np.unique(predictions):
    #        df[label] = np.nan
    #return df
    df     = df[["Unrelated", "Third", "Second", "First", "Self"]]
    return df.reindex(["Unrelated", "Third", "Second", "First", "Self"])


def main():

    (truths, preds) = ([], [])
    for line in sys.stdin:
        line = line.strip('\n').split(" ")
        truths.append(line[0])
        preds.append(grups_to_degree(line[1]))

    #print_confusion_matrix(make_confusion_matrix(truths, preds))
    print(pd_print_confusion_matrix(truths, preds))

if __name__ == '__main__':
    main()
