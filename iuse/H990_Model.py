import argparse
import sys, os, platform
current_directory = os.path.dirname(os.path.abspath(__file__))
root_path = os.path.abspath(os.path.dirname(current_directory) + os.path.sep + ".")
sys.path.append(root_path)
import pickle
from descnucleotide import PSTNPss,binary,Pse
from pubscripts import read_fasta_sequences, save_file
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from featureselection import select_features
import warnings
warnings.filterwarnings("ignore")
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
data_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\data' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/data'
data_path1 = os.path.abspath(
        os.path.dirname(pPath) + os.path.sep + ".") + r'\iuse' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/iuse'

def Standard_feature(X_st):
    scaler =MinMaxScaler()
    X_st = scaler.fit_transform(X_st)
    return X_st


def read_feature_selection(file):
    result = []

    with open(data_path1 + '/'+file, 'r')as f:
        for line in f.readlines():
            result.append(list(map(str, line.strip().split(','))))
    return result
didna_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content', 'A-philicity',
              'Propeller twist', 'Duplex stability:(freeenergy)',
              'Duplex tability(disruptenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
              'Stabilising energy of Z-DNA', 'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH',
              'Breslauer_dS', 'Electron_interaction', 'Hartman_trans_free_energy', 'Helix-Coil_transition',
              'Ivanov_BA_transition', 'Lisser_BZ_transition', 'Polar_interaction', 'SantaLucia_dG',
              'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Stacking_energy',
              'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS', 'Watson-Crick_interaction', 'Twist', 'Tilt', 'Roll',
              'Shift', 'Slide', 'Rise',
              'Clash Strength', 'Roll_roll', 'Twist stiffness', 'Tilt stiffness', 'Shift_rise',
              'Adenine content', 'Direction', 'Twist_shift', 'Enthalpy1', 'Twist_twist', 'Roll_shift',
              'Shift_slide', 'Shift2', 'Tilt3', 'Tilt1', 'Tilt4', 'Tilt2', 'Slide (DNA-protein complex)1',
              'Tilt_shift', 'Twist_tilt', 'Twist (DNA-protein complex)1', 'Tilt_rise', 'Roll_rise',
              'Stacking energy', 'Stacking energy1', 'Stacking energy2', 'Stacking energy3', 'Propeller Twist',
              'Roll11', 'Rise (DNA-protein complex)', 'Tilt_tilt', 'Roll4', 'Roll2', 'Roll3', 'Roll1',
              'Minor Groove Size', 'GC content', 'Slide_slide', 'Enthalpy', 'Shift_shift', 'Slide stiffness',
              'Melting Temperature1', 'Flexibility_slide', 'Minor Groove Distance',
              'Rise (DNA-protein complex)1', 'Tilt (DNA-protein complex)', 'Guanine content',
              'Roll (DNA-protein complex)1', 'Entropy', 'Cytosine content', 'Major Groove Size', 'Twist_rise',
              'Major Groove Distance', 'Twist (DNA-protein complex)', 'Purine (AG) content',
              'Melting Temperature', 'Free energy', 'Tilt_slide', 'Major Groove Width', 'Major Groove Depth',
              'Wedge', 'Free energy8', 'Free energy6', 'Free energy7', 'Free energy4', 'Free energy5',
              'Free energy2', 'Free energy3', 'Free energy1', 'Twist_roll', 'Shift (DNA-protein complex)',
              'Rise_rise', 'Flexibility_shift', 'Shift (DNA-protein complex)1', 'Thymine content', 'Slide_rise',
              'Tilt_roll', 'Tip', 'Keto (GT) content', 'Roll stiffness', 'Minor Groove Width', 'Inclination',
              'Entropy1', 'Roll_slide', 'Slide (DNA-protein complex)', 'Twist1', 'Twist3', 'Twist2', 'Twist5',
              'Twist4', 'Twist7', 'Twist6', 'Tilt (DNA-protein complex)1', 'Twist_slide', 'Minor Groove Depth',
              'Roll (DNA-protein complex)', 'Rise2', 'Persistance Length', 'Rise3', 'Shift stiffness',
              'Probability contacting nucleosome core', 'Mobility to bend towards major groove', 'Slide3',
              'Slide2', 'Slide1', 'Shift1', 'Bend', 'Rise1', 'Rise stiffness',
              'Mobility to bend towards minor groove']

tridna_list = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
               'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
               'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

dirna_list = ['Slide (RNA)', 'Adenine content', 'Hydrophilicity (RNA)', 'Tilt (RNA)', 'Stacking energy (RNA)',
              'Twist (RNA)', 'Entropy (RNA)', 'Roll (RNA)', 'Purine (AG) content', 'Hydrophilicity (RNA)1',
              'Enthalpy (RNA)1', 'GC content', 'Entropy (RNA)1', 'Rise (RNA)', 'Free energy (RNA)',
              'Keto (GT) content', 'Free energy (RNA)1', 'Enthalpy (RNA)', 'Guanine content', 'Shift (RNA)',
              'Cytosine content', 'Thymine content']

myDict = {
    'DAC': {'DNA': didna_list, 'RNA': dirna_list},
    'DCC': {'DNA': didna_list, 'RNA': dirna_list},
    'DACC': {'DNA': didna_list, 'RNA': dirna_list},
    'TAC': {'DNA': tridna_list, 'RNA': []},
    'TCC': {'DNA': tridna_list, 'RNA': []},
    'TACC': {'DNA': tridna_list, 'RNA': []},
    'PseDNC': {'DNA': didna_list, 'RNA': dirna_list},
    'PseKNC': {'DNA': didna_list, 'RNA': dirna_list},
    'PCPseDNC': {'DNA': didna_list, 'RNA': dirna_list},
    'PCPseTNC': {'DNA': tridna_list, 'RNA': []},
    'SCPseDNC': {'DNA': didna_list, 'RNA': dirna_list},
    'SCPseTNC': {'DNA': tridna_list, 'RNA': []},
}

myDictDefault = {
    'DAC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
            'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'DCC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
            'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'DACC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
             'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'TAC': {'DNA': ['Dnase I', 'Bendability (DNAse)'],
            'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'TCC': {'DNA': ['Dnase I', 'Bendability (DNAse)'],
            'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'TACC': {'DNA': ['Dnase I', 'Bendability (DNAse)'],
             'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'PseDNC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
               'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'PseKNC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
               'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'PCPseDNC': {
        'DNA': ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'A-philicity', 'Propeller twist',
                'Duplex stability:(freeenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
                'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH', 'Electron_interaction',
                'Hartman_trans_free_energy', 'Helix-Coil_transition', 'Lisser_BZ_transition', 'Polar_interaction',
                'SantaLucia_dG', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Sugimoto_dG', 'Sugimoto_dH',
                'Sugimoto_dS', 'Duplex tability(disruptenergy)', 'Stabilising energy of Z-DNA', 'Breslauer_dS',
                'Ivanov_BA_transition', 'SantaLucia_dH', 'Stacking_energy', 'Watson-Crick_interaction',
                'Dinucleotide GC Content', 'Twist', 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise'],
        'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'PCPseTNC': {'DNA': ['Dnase I', 'Bendability (DNAse)'], 'RNA': []},
    'SCPseDNC': {'DNA': ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist'],
                 'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},
    'SCPseTNC': {'DNA': ['Dnase I', 'Bendability (DNAse)'], 'RNA': []},
}

myKmer = {
    'DAC': 2, 'DCC': 2, 'DACC': 2,
    'TAC': 3, 'TCC': 3, 'TACC': 3
}

myDataFile = {
    'DAC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'DCC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'DACC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'TAC': {'DNA': 'tridnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'TCC': {'DNA': 'tridnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'TACC': {'DNA': 'tridnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'PseDNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'PseKNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'PCPseDNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'PCPseTNC': {'DNA': 'tridnaPhyche.data', 'RNA': ''},
    'SCPseDNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},
    'SCPseTNC': {'DNA': 'tridnaPhyche.data', 'RNA': ''},
}

def check_Pse_arguments(method,type):

    myNum = 0

    myIndex = []
    myProperty = {}
    dataFile = ''

    if myNum == 0:
        myIndex = myDictDefault[method][type]
        dataFile = myDataFile[method][type]
        # print(myIndex,dataFile)
    if dataFile != '':
        with open(data_path + '/' + dataFile, 'rb') as f:
        # with open(dataFile, 'rb') as f:
            myProperty = pickle.load(f)

    if len(myIndex) == 0 or len(myProperty) == 0:
        print('Error: arguments is incorrect.')
        sys.exit(1)

    return myIndex, myProperty#, args.lamadaValue, args.weight, args.kmer


def ind_fs(X_,en,pe):
    if en == 'binary':
        feature_num7 = pe
        feature_selection7 = read_feature_selection('feature_selection/feature_selection7.txt')
        X_binary = binary.binary(X_)
        X_binary_fs = select_features.select_features(X_binary,feature_selection7, feature_num7)
        X_binary_fs = np.array(X_binary_fs[1:])
        return X_binary_fs
    elif en == 'PSTNPss':
        feature_num13 = pe  # int((len(Bm_selection13[1]) - 2) * pe)
        feature_selection13 = read_feature_selection('feature_selection/feature_selection13.txt')
        X_pstnpss = PSTNPss.PSTNPss(X_,'H')
        X_pstnpss_fs = select_features.select_features(X_pstnpss,feature_selection13, feature_num13)
        X_pstnpss_fs = np.array(X_pstnpss_fs[1:])
        return X_pstnpss_fs
    elif en== 'PseKNC':
        my_property_name, my_property_value = check_Pse_arguments(method='PseKNC', type='RNA')
        feature_num15 = pe
        feature_selection15 = read_feature_selection('feature_selection/feature_selection15.txt')
        X_PSEKNC = Pse.make_PseKNC_vector(X_, my_property_name, my_property_value, 2, 0.1, 3)
        X_PSEKNC_FS = select_features.select_features(X_PSEKNC,feature_selection15, feature_num15)
        X_PSEKNC_FS = np.array(X_PSEKNC_FS[1:])
        return X_PSEKNC_FS

def predictor(inputfile, outputfile):
    fastas = read_fasta_sequences.read_nucleotide_sequences(inputfile)
    fastas=np.array(fastas)
    X = fastas[:, 0:]
    X_per1 = ind_fs(X, 'PSTNPss', 19)
    X_per2 = ind_fs(X, 'PSTNPss', 19)
    X_per3 = ind_fs(X, 'PseKNC', 64)
    X_per4 = ind_fs(X, 'binary', 7)
    X_per = np.array(np.concatenate(tuple((X_per1, X_per2[:, 2:], X_per3[:, 2:], X_per4[:, 2:])), axis=1))
    X_per = Standard_feature(X_per[:, 2:])
    modelfile='Pickle/H_model.pickle'
    with open(data_path1+'/'+modelfile, 'rb') as f:
        sclf_restore = pickle.load(f)
    scores = sclf_restore.predict_proba(X_per)
    #
    prdiction_result_with_suquence=[]
    pre_result=np.zeros((len(X),4)).astype(np.str)
    pre_result[:,0],pre_result[:,1],pre_result[:,2:]=fastas[:,0],fastas[:,1],scores
    prdiction_result_with_suquence.append(pre_result)
    save_file.save_CV_result_binary1(prdiction_result_with_suquence,outputfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="input fasta file")
    parser.add_argument("--o", required=True, help="output file")
    args = parser.parse_args()
    predictor(args.i,args.o)
