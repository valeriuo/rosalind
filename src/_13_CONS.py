# Consensus and Profile
# http://rosalind.info/problems/cons/
import sys


def main():
    profile = {'A':[], 'C':[], 'G':[], 'T':[], 'N':[]}
    consensus = ''
    base_count = 0

    with open('../data/rosalind_cons.txt', 'r') as fp:
        seq = ''
        for line in fp:
            if line[0] != '>':
                seq += line.replace('\n', '')
            else:
                if base_count == 0 and len(seq) > 0:
                    base_count = len(seq)

                    profile['A'] = [0]*base_count
                    profile['C'] = [0]*base_count
                    profile['G'] = [0]*base_count
                    profile['T'] = [0]*base_count
                    profile['N'] = [0]*base_count                    
                    
                for i, b in enumerate(seq):
                    l = profile.get(b, profile['N'])
                    if i < base_count:
                        l[i] += 1       

                seq = ''
        # process last sequence        
        for i, b in enumerate(seq):
            l = profile.get(b, profile['N'])
            if i < base_count:
                l[i] += 1  
                
    lA = profile['A']
    lC = profile['C']
    lG = profile['G']
    lT = profile['T']
    
    for i in range(base_count):
        li = [(lA[i],'A'), (lC[i],'C'), (lT[i],'G'), (lT[i],'T')]
        lb = [lA[i], lC[i], lG[i], lT[i]]
        consensus = consensus + li[lb.index(max(lb))][1]
        
    print(consensus)
    for i in profile.keys():
        if i != 'N':
            print(i, end=': ')
            for j in profile[i]:
                print(j, end=' ')
            print('')
    
if __name__ == '__main__':
    main()