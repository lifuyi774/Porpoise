import argparse
import sys, os, platform
current_directory = os.path.dirname(os.path.abspath(__file__))
root_path = os.path.abspath(os.path.dirname(current_directory) + os.path.sep + ".")
sys.path.append(root_path)
import pickle
from descnucleotide import PSTNPss,NCP
from pubscripts import read_fasta_sequences,save_file
import numpy as np
from sklearn.preprocessing import Normalizer
from featureselection import select_features
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
data_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\data' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/data'
data_path1 = os.path.abspath(
        os.path.dirname(pPath) + os.path.sep + ".") + r'\iuse' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/iuse'

def Standard_feature(X_st):
    scaler =Normalizer()
    X_st= scaler.fit_transform(X_st)
    return X_st

def read_feature_selection(file):
    result = []
    # with open(file, 'r')as f:
    with open(data_path1 + '/' + file, 'r')as f:
        for line in f.readlines():
            result.append(list(map(str, line.strip().split(','))))
    return result

def train_valid_fs(X_,en,pe):
    if en =='NCP':
        feature_num9 = pe
        feature_selection9 = read_feature_selection('feature_selection/feature_selection9-2.txt')
        X_ncp = NCP.NCP(X_)
        X_ncp_fs = select_features.select_features(X_ncp, feature_selection9, feature_num9)
        X_ncp_fs = np.array(X_ncp_fs[1:])
        return X_ncp_fs
    elif en =='PSTNPss':
        feature_num13 = pe
        feature_selection13 = read_feature_selection('feature_selection/feature_selection13-2.txt')
        X_pstnpss = PSTNPss.PSTNPss(X_,'M')
        X_pstnpss_fs = select_features.select_features(X_pstnpss,  feature_selection13, feature_num13)
        X_pstnpss_fs = np.array(X_pstnpss_fs[1:])
        return X_pstnpss_fs


def predictor(inputfile, outputfile):
    fastas = read_fasta_sequences.read_nucleotide_sequences(inputfile)
    fastas = np.array(fastas)
    X = fastas[:, 0:]
    X_per1 = train_valid_fs(X,'PSTNPss', 15)
    X_per2 = train_valid_fs(X, 'NCP', 35)
    X_per1 = Standard_feature(X_per1[:, 2:])
    X_per2 = Standard_feature(X_per2[:, 2:])
    X_per = np.array(np.concatenate(tuple((X_per1, X_per2)), axis=1))

    modelfile = 'Pickle/M_model.pickle'
    with open(data_path1 + '/' + modelfile, 'rb') as f:
        sclf_restore = pickle.load(f)
    scores = sclf_restore.predict_proba(X_per)
    prdiction_result_with_suquence = []
    pre_result = np.zeros((len(X), 4)).astype(np.str)
    pre_result[:, 0], pre_result[:, 1], pre_result[:, 2:] = fastas[:, 0], fastas[:, 1], scores
    prdiction_result_with_suquence.append(pre_result)
    save_file.save_CV_result_binary1(prdiction_result_with_suquence, outputfile)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="input fasta file")
    parser.add_argument("--o", required=True, help="output file")
    args = parser.parse_args()
    predictor(args.i,args.o)

