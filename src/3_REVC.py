# DNA reverse complement generator
import sys

def revc(input):
    output = ""
    for c in input:
        if c == 'A':
            output = "T" + output
        elif c == 'C':
            output = "G" + output
        elif c == 'G':
            output = "C" + output
        elif c == 'T':
            output = "A" + output
            
    return output 

def main():
    with open('../data/rosalind_revc.txt', 'r') as myfile:
        source=myfile.read().replace('\n', '')
    print(revc(source))

if __name__ == '__main__':
    main()