# DNA base counter
import sys

def dna(input):
    nA = nC = nG = nT = nO = 0
    for c in input:
        if c == 'A':
            nA += 1
        elif c == 'C':
            nC += 1
        elif c == 'G':
            nG += 1
        elif c == 'T':
            nT += 1
        else:
            nO += 1
    output = "{} {} {} {}".format(nA, nC, nG, nT)
    
    return output    

def main():
    with open(sys.argv[1], 'r') as myfile:
        source=myfile.read().replace('\n', '')
    print(dna(source))
    
if __name__ == '__main__':
    main()