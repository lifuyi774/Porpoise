import re, os, sys
import argparse
def processing(file,outfile,nb_windows):
    fw = open(outfile,'w')
    if os.path.exists(file) == False:
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACGTU-]', '-', ''.join(array[1:]).upper())
        sequence = re.sub('U', 'T', sequence)
        nb_windows = int(nb_windows)
        empty_aa = '-'
        focus =['T']
        for pos in range(len(sequence)):
            mid_aa = sequence[pos]
            if not (mid_aa in focus):
                continue

            start = 0
            if pos - nb_windows + 1 > 0:
                start = pos - nb_windows
            left_seq = sequence[start:pos]

            end = len(sequence)
            if pos + nb_windows < end:
                end = pos + nb_windows + 1
            right_seq = sequence[pos + 1:end]

            if len(left_seq) < nb_windows:
                if empty_aa is None:
                    continue
                nb_lack = nb_windows - len(left_seq)
                left_seq = ''.join([empty_aa for _count in range(nb_lack)]) + left_seq

            if len(right_seq) < nb_windows:
                if empty_aa is None:
                    continue
                nb_lack = nb_windows - len(right_seq)
                right_seq = right_seq + ''.join([empty_aa for _count in range(nb_lack)])
            final_seq = left_seq + mid_aa + right_seq
            fw.write('>' + header + '\t'+str(pos+1) + '\n' + final_seq + '\n')
    fw.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="input fasta file")
    parser.add_argument("--o", required=True, help="file of processed data")
    parser.add_argument("--w", required=True, help="sequence windows(10 or 15)")
    args = parser.parse_args()
    processing(args.i,args.o,args.w)