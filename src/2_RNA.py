# RNA transcript of a DNA sequence
import sys

def rna(input):
    output = ""
    for c in input:
        if c == 'T':
            output += "U"
        else:
            output += c
            
    return output 

def main():
    with open('../data/rosalind_rna.txt', 'r') as myfile:
        source=myfile.read().replace('\n', '')
    print(rna(source))

if __name__ == '__main__':
    main()