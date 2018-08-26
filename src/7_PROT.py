# RNA to protein translation
import sys

def amino_acid(input):
    rna_codon_table = { 'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V', 'UUC':'F', 'CUC':'L', 'AUC':'I', 'GUC':'V', 'UUA':'L', 'CUA':'L', 'AUA':'I', 'GUA':'V', 'UUG':'L', 'CUG':'L', 'AUG':'M', 'GUG':'V', 'UCU':'S', 'CCU':'P', 'ACU':'T', 'GCU':'A', 'UCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A', 'UCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A', 'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A', 'UAU':'Y', 'CAU':'H', 'AAU':'N', 'GAU':'D', 'UAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D', 'UAA':'Stop', 'CAA':'Q', 'AAA':'K', 'GAA':'E', 'UAG':'Stop', 'CAG':'Q', 'AAG':'K', 'GAG':'E', 'UGU':'C', 'CGU':'R', 'AGU':'S', 'GGU':'G', 'UGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G', 'UGA':'Stop', 'CGA':'R', 'AGA':'R', 'GGA':'G', 'UGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G' }
    
    return rna_codon_table.get(input)

def main():
    with open('../data/rosalind_prot.txt', 'r') as myfile:
        rna=myfile.read().replace('\n', '')
    #codons = [rna[3*i:3*(i+1)] for i in range(len(rna)//3)]
    amino_acids = [amino_acid(rna[3*i:3*(i+1)]) for i in range(len(rna)//3) if amino_acid(rna[3*i:3*(i+1)]) is not None]
    #print(codons)
    #print(amino_acids)
    protein = ''.join(amino_acids[:amino_acids.index('Stop')])
    print(protein)

if __name__ == '__main__':
    main()