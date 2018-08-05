# Finding a Shared Motif
# http://rosalind.info/problems/lcsm/
import sys

def lcsm(r1, r2):
    cs = []
    for i in range(len(r1)-1):
        for j in range(2, len(r1)-i+1):
            if (r1[i:i+j] in r2):
                if (r1[i:i+j] not in cs):
                    cs.append(r1[i:i+j])
                if (j > 2):
                    cs.remove(r1[i:i+j-1])
            else:
                break
    return cs

def main():

    reads = []
    lcs = []
    with open('../data/rosalind_lcsm.txt', 'r') as fp:
        seq = ''
        for line in fp:
            if line[0] != '>':
                seq += line.replace('\n', '')
            else:
                if len(seq) > 0:
                    reads.append(seq)
                    seq = ''

        # process last sequence        
        reads.append(seq)
        
    lcs = lcs + lcsm(reads[0], reads[1])
    lcs.sort(key = len, reverse=True)
    
    for i in range(2, len(reads)):
        for l in lcs:
            if l not in reads[i]:
                lcs.remove(l)
            else:
                break
                
    print(lcs[0])
                    
    
if __name__ == '__main__':
    main()