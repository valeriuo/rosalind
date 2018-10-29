#Enumerating k-mers Lexicographically
import sys
    
def main():
    result = ['']
    with open('../data/rosalind_lexf.txt', 'r') as fp:
        alphabet = fp.readline().split()
        n = int(fp.readline())
        
        alphabet.sort()
        for i in range(n):
            result = [s1+s2 for s1 in result for s2 in alphabet]
    print(*result, sep='\n')
    
if __name__ == '__main__':
    main()