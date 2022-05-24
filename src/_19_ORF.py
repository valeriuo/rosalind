# Open Reading Frames
import sys
from _3_REVC import revc

def dna2aa(input):
    dna_codon_table = { 'TTT':'F', 'CTT':'L', 'ATT':'I', 'GTT':'V', 'TTC':'F', 'CTC':'L', 'ATC':'I', 'GTC':'V', 'TTA':'L', 'CTA':'L', 'ATA':'I', 'GTA':'V', 'TTG':'L', 'CTG':'L', 'ATG':'M', 'GTG':'V', 'TCT':'S', 'CCT':'P', 'ACT':'T', 'GCT':'A', 'TCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A', 'TCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A', 'TCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A', 'TAT':'Y', 'CAT':'H', 'AAT':'N', 'GAT':'D', 'TAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D', 'TAA':'Stop', 'CAA':'Q', 'AAA':'K', 'GAA':'E', 'TAG':'Stop', 'CAG':'Q', 'AAG':'K', 'GAG':'E', 'TGT':'C', 'CGT':'R', 'AGT':'S', 'GGT':'G', 'TGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G', 'TGA':'Stop', 'CGA':'R', 'AGA':'R', 'GGA':'G', 'TGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G' }
    
    return dna_codon_table.get(input)
    
def peptides(input):
    output = set()
    amino_acids = [dna2aa(input[3*i:3*(i+1)]) for i in range(len(input)//3) if dna2aa(input[3*i:3*(i+1)]) is not None]
    #print(input)
    #print(amino_acids)

    for i in range(len(amino_acids)):
        if amino_acids[i] == 'M' and 'Stop' in amino_acids[i+1:]:
            output.add(''.join(amino_acids[i:amino_acids.index('Stop',i+1)]))
            
    return output

def main():
    if len(sys.argv) < 2:
        data_file = '../data/rosalind_orf.txt'
    else:
        data_file = sys.argv[1]

    with open(data_file, 'r') as fp:
        fseq = ''
        for line in fp:
            if line[0] != '>':
                fseq += line.replace('\n', '')
            else:
                if len(fseq) > 0:
                    # process sequence
                    all = peptides(fseq)
                    all |= peptides(fseq[1:])
                    all |= peptides(fseq[2:])
                    
                    rseq = revc(fseq)
                    all |= peptides(rseq)
                    all |= peptides(rseq[1:])
                    all |= peptides(rseq[2:])
                    
                    print('\n'.join(all))
                    # clear sequence
                    fseq = ''

        # process last sequence           
        #fseq = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
        all = peptides(fseq)
        all |= peptides(fseq[1:])
        all |= peptides(fseq[2:])
        
        rseq = revc(fseq)
        all |= peptides(rseq)
        all |= peptides(rseq[1:])
        all |= peptides(rseq[2:])
        
        print('\n'.join(all))

if __name__ == '__main__':
    main()