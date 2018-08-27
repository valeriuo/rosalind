# RNA Splicing
import sys
import re
from _19_ORF import dna2aa

def splice(input, intron):
    output = ''
    intron_pos = [(m.start(), m.end()) for m in re.finditer(intron, input)]
    start = 0
    for i in intron_pos:
        output += input[start:i[0]]
        start = i[1]
    output += input[start:]
        
    return output


def main():
    with open('../data/rosalind_splc.txt', 'r') as fp:
        seq = ''
        read_seq = 0
        intron = ''
        for line in fp:
            if line[0] != '>':
                if read_seq is 0:
                    seq += line.replace('\n', '')
                else:
                    intron += line.replace('\n', '')
            else:
                if len(seq) > 0:
                    if read_seq is 0:
                        read_seq = 1
                    # process sequence
                    if len(intron) > 0:
                        seq = splice(seq, intron)
                        # clear sequence
                        intron = ''

        # process last sequence           
        if len(seq) > 0 and len(intron) > 0:
            seq = splice(seq, intron)
            
        amino_acids = [dna2aa(seq[3*i:3*(i+1)]) for i in range(len(seq)//3) if dna2aa(seq[3*i:3*(i+1)]) is not None]
        protein = ''.join(amino_acids[amino_acids.index('M'):amino_acids.index('Stop')])
        print(protein)

if __name__ == '__main__':
    main()