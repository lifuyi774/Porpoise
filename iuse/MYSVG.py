#!/usr/bin/env python
# _*_coding:utf-8_*_
import argparse
header = """<?xml version="1.0" standalone="no"?>

<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

<svg width="100%" height="100%" version="1.1" xmlns="http://www.w3.org/2000/svg">
"""

tail = """
</svg>
"""


def dropProtein(proteinName, sequenceLength, Y):
    rect = '<rect x="200" y="%d" rx="5" ry="5" width="800" height="10" style="fill:#CCCCFF;stroke:black; stroke-width:1;opacity:0.5"/>' % Y
    textProtein = '<text x="50" y="%d" style="font-weight:bold;font-size:14pt;" fill="#000000">%s</text>' % (
    Y + 10, proteinName)
    beginLine = '<line x1="200" y1="%d" x2="200" y2="%d" style="stroke:#000066;stroke-width:1"/>' % (Y - 10, Y + 20)
    endLine = '<line x1="1000" y1="%d" x2="1000" y2="%d" style="stroke:#000066;stroke-width:1"/>' % (Y - 10, Y + 20)
    beginText = '<text x="180" y="%d" fill="#000000">1</text>' % (Y + 10)
    endText = '<text x="1010" y="%d" fill="#000000">%d</text>' % (Y + 10, sequenceLength)
    return rect + '\n' + textProtein + '\n' + beginLine + '\n' + endLine + '\n' + beginText + '\n' + endText + '\n'


def labelUp(position, sequenceLength, Y, X=200):
    x = position / sequenceLength * 800 + X
    polyline = '<polyline points="%d,%d %d,%d %d,%d" style="fill:white;stroke:#3333CC;stroke-width:1"/>' % (
        x, Y, x, Y - 15, x + 30, Y - 25)
    circle = '<circle cx="%d" cy="%d" r="5" stroke="#FF0033" stroke-width="1" fill="#FFFFFF"/>' % (x + 30, Y - 25)
    text = '<text x="%d" y="%d" fill="#3333FF">U%d</text>' % (x + 20, Y - 35, position)
    return polyline + '\n' + circle + '\n' + text


def labelDown(position, sequenceLength, Y, X=200):
    x = position / sequenceLength * 800 + X
    Y = Y + 10
    polyline = '<polyline points="%d,%d %d,%d %d,%d" style="fill:white;stroke:#3333CC;stroke-width:1"/>' % (
        x, Y, x, Y + 15, x + 30, Y + 25)
    circle = '<circle cx="%d" cy="%d" r="5" stroke="#FF0033" stroke-width="1" fill="#FFFFFF"/>' % (x + 30, Y + 25)
    text = '<text x="%d" y="%d" fill="#3333FF">U%d</text>' % (x + 20, Y + 45, position)
    return polyline + '\n' + circle + '\n' + text


def createBarChart(rawdatafile, threshold, modelout, species, output):
    modelThreshold = {
        'H.sapiens': [0.221,0.267,0.781],
        'S.cerevisiae': [0.066,0.068,0.914],
        'M.musculus': [0.256,0.315,0.788]
    }
    myDict = { 'LOW': 0, 'MEDIUM': 1, 'HIGH': 2}


    rnaLength = {}
    with open(rawdatafile) as f:
        mySequence = f.read()

    array = mySequence.split('>')[1:]
    for i in array:
        fasta = i.strip().replace(' ', '').split('\n')
        rnaLength[fasta[0]] = len(''.join(fasta[1:]))

    na_pos_sc = {}
    with open(modelout) as f1:
        name_sequence_score = f1.read()
    result_arr = name_sequence_score.split('\n')
    result_arr.remove('')

    for result in result_arr:
        re_arr = result.split('\t')
        if len(re_arr) != 3:
            pass
        name_pos, score = re_arr[0], float(re_arr[2])
        name_pos_arr = name_pos.split('_')
        name = name_pos_arr[0]

        if name not in na_pos_sc.keys():
            na_pos_sc[name] = {}
            na_pos_sc[name][name_pos_arr[1]]=score
        else:
            na_pos_sc[name][name_pos_arr[1]]=score


    rna = {}
    rnaOrder = []
    for na in na_pos_sc.keys():
        if na in rnaOrder:
            pass
        else:
            rnaOrder.append(na)
            rna[na] = []
        for p in na_pos_sc[na].keys():
            if na_pos_sc[na][p] >= modelThreshold[species][myDict[threshold]]:
                rna[na].append(int(p))


    with open(output, 'w') as f:
        f.write(header)
        f.write(
            '<text x="400" y="50" style="font-weight:bold;font-size:20pt;" fill="#3366CC">Distribution of Pseudouridine sites</text>')
        f.write(
            '<text x="450" y="90" style="font-weight:bold;font-size:16pt;" fill="#808080">Threshold: %s</text>' % threshold)
        y = 200
        for p in rnaOrder:
            f.write(dropProtein(p, rnaLength[p], y))
            tag = 0
            for pos in rna[p]:
                if tag % 2 == 0:
                    f.write(labelUp(pos, rnaLength[p], y))
                else:
                    f.write(labelDown(pos, rnaLength[p], y))
                tag = tag + 1
            y = y + 150

        f.write(tail)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="raw data file")
    parser.add_argument("--t", required=True, help="threshold(LOW,MEDIUM,HIGH)")
    parser.add_argument("--m", required=True, help="model's result file")
    parser.add_argument("--s", required=True, help="species(H.sapiens,S.cerevisiae,M.musculus)")
    parser.add_argument("--o", default='result.svg', help="result file ")
    args = parser.parse_args()
    createBarChart(args.i, args.t, args.m, args.s, args.o)




