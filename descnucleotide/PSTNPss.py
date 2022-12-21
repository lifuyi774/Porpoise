#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys, os, platform
import re
import numpy as np
from pubscripts import *
import csv
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import check_sequences


def PSTNPss(fastas,typ, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "PSTNP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    for i in fastas:
        if re.search('[^ACGT-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT[U]" are allowed by this encoding scheme.')
            return 0

    encodings = []
    header = ['#', 'label']
    for pos in range(len(fastas[0][1])-2):
        header.append('Pos.%d' %(pos+1))
    encodings.append(header)

    nucleotides = ['A', 'C', 'G', 'T']
    trinucleotides = [n1 + n2 + n3 for n1 in nucleotides for n2 in nucleotides for n3 in nucleotides]
    order = {}
    for i in range(len(trinucleotides)):
        order[trinucleotides[i]] = i
    if typ == 'H':
        matrix_po = np.fromfile(father_path+'/'+'H_po.bin', dtype=float)
        matrix_po = matrix_po.reshape(19, 64)
        matrix_ne = np.fromfile(father_path + '/' + 'H_ne.bin', dtype=float)
        matrix_ne = matrix_ne.reshape(19, 64)
        positive_number = 4327
        negative_number = 3732

    elif typ == 'S':
        matrix_po = np.fromfile(father_path + '/' + 'S_po.bin', dtype=float)
        matrix_po = matrix_po.reshape(29, 64)
        matrix_ne = np.fromfile(father_path + '/' + 'S_ne.bin', dtype=float)
        matrix_ne = matrix_ne.reshape(29, 64)
        positive_number = 1147
        negative_number = 1147

    elif typ == 'M':
        matrix_po = np.fromfile(father_path + '/' + 'M_po.bin', dtype=float)
        matrix_po = matrix_po.reshape(19, 64)
        matrix_ne = np.fromfile(father_path + '/' + 'M_ne.bin', dtype=float)
        matrix_ne = matrix_ne.reshape(19, 64)
        positive_number = 3174
        negative_number = 3174

    for i in fastas:

        name, sequence= i[0], i[1]
        label='0'
        code = [name, label]

        for j in range(len(sequence) - 2):
            if re.search('-', sequence[j: j+3]):
                code.append(0)
            else:
                p_num, n_num = positive_number, negative_number
                po_number = matrix_po[j][order[sequence[j: j+3]]]
                ne_number = matrix_ne[j][order[sequence[j: j+3]]]
                code.append(po_number/p_num - ne_number/n_num)
        encodings.append(code)
    return encodings


