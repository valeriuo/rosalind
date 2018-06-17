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

def rna(input):
    output = ""
    for c in input:
        if c == 'T':
            output += "U"
        else:
            output += c
            
    return output 
    
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
    with open(sys.argv[1], 'r') as myfile:
        source=myfile.read().replace('\n', '')
    print(revc(source))

if __name__ == '__main__':
    main()