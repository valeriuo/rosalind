# Locating Restriction Sites
import sys
from _3_REVC import revc

def is_revp(input):
    l = len(input)

    return input[:l//2] == revc(input[(l+1)//2:])

    
def gen_revp(input):
    l = len(input)
    for i in range(l-3):
        for j in range(4, min(12, l-i) + 1, 2):
            if is_revp(input[i:i+j]):
                yield (i+1, j)

def main():
    with open('../data/rosalind_revp.txt', 'r') as fp:
        seq = ''
        for line in fp:
            if line[0] != '>':
                seq += line.replace('\n', '')
            else:
                if len(seq) > 0:
                    # process sequence
                    for r in gen_revp(seq):
                        print(*r)
                    # clear sequence
                    seq = ''

        # process last sequence
        for r in gen_revp(seq):
            print(*r)

if __name__ == '__main__':
    main()