# Transition/transversion ration between two strings of equal length
import sys

def tran_ratio(s, t):
    trs = 0
    trv = 0

    for cs, ct in zip(s,t):
        if cs != ct:
            if set([cs, ct]) == set(['A', 'G']) or set([cs, ct]) == set(['C', 'T']):
                trs += 1
            else:
                trv += 1

    if trv > 0:
        return trs/trv

    return trs        
    
def main():
    if len(sys.argv) < 2:
        data_file = '../data/rosalind_tran.txt'
    else:
        data_file = sys.argv[1]

    s1 = ''
    s2 = ''
    with open(data_file, 'r') as fp:
        fseq = ''
        for line in fp:
            if line[0] != '>':
                fseq += line.replace('\n', '')
            else:
                s1 = fseq
                fseq = ''
        s2 = fseq                

    if len(s1) != len(s2):
        print (f"[error]: Strings {s1} and {s2} don't have the same length")
    else:
        print(tran_ratio(s1,s2))
if __name__ == '__main__':
    main()