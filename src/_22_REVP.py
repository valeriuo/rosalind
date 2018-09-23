# Locating Restriction Sites
import sys
from _3_REVC import revc

def is_revp(input):
    l = len(input)
    
    if input[:l//2] == revc(input[(l+1)//2:]):
        return True
    return False

def main():
    with open('../data/rosalind_revp.txt', 'r') as fp:
        seq = ''
        for line in fp:
            if line[0] != '>':
                seq += line.replace('\n', '')
            else:
                if len(seq) > 0:
                    # process sequence
                    l = len(seq)
                    for i in range(l-3):
                        for j in range(4, min(12, l-i) + 1, 2):
                            if is_revp(seq[i:i+j]):
                                print(i+1, j)
                    # clear sequence
                    seq = ''

        # process last sequence
        l = len(seq)
        for i in range(l-3):
            for j in range(4, min(12, l-i) + 1, 2):
                if is_revp(seq[i:i+j]):
                    print(i+1, j)

if __name__ == '__main__':
    main()