# Overlap Graphs
# http://rosalind.info/problems/grph/
import sys

def is_connected(r1, r2, k):
    if (r1[-k:] == r2[:k]):
        return 1
    elif (r2[-k:] == r1[:k]):
        return -1
    else: 
        return 0


def main():

    names = []
    reads = []
    with open('../data/rosalind_grph.txt', 'r') as fp:
        seq = ''
        for line in fp:
            if line[0] != '>':
                seq += line.replace('\n', '')
            else:
                names.append(line[1:].replace('\n', ''))
                if len(seq) > 0:
                    reads.append(seq)
                    seq = ''

        # process last sequence        
        reads.append(seq)
        seq = ''
        
    for i in range(len(reads)):
        for j in range(i+1, len(reads)):
            if is_connected(reads[i], reads[j], 3) == 1:
                print(names[i], names[j])
            elif is_connected(reads[i], reads[j], 3) == -1:
                print(names[j], names[i])
    
if __name__ == '__main__':
    main()