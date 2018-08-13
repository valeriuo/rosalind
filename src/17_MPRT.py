# Finding a Protein Motif
# http://rosalind.info/problems/mprt/
import sys
import urllib.request as ur
import re

def mprt(name, seq):
    pat = re.compile('(?=(N[^P][ST][^P]))') # use a lookahead capturing group to catch the overlapping motifs
    coord = [m.span()[0]+1 for m in list(pat.finditer(seq))]
    if len(coord) > 0:
        print(name)
        print(*coord, sep=' ')
    
def main():

    with open('../data/rosalind_mprt.txt', 'r') as prot_names:
        for line in prot_names:
            name = line.replace('\n', '')    
            with ur.urlopen('https://www.uniprot.org/uniprot/' + name + '.fasta') as fp:
                seq = ''
                for line in fp:
                    line = line.decode('utf-8')
                    if line[0] != '>':
                        seq += line.replace('\n', '')
                    else:
                        if len(seq) > 0:
                            mprt(name, seq)
                            seq = ''
            # process last sequence        
            mprt(name, seq)
    
if __name__ == '__main__':
    main()