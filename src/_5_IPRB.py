# Mendel's First Law
import sys

def iprb(k, m, n):
    if k + m + n < 2 :
        return 0
    return (k*(k-1)/2 + k*m + m*(m-1)*3/8 + m*n/2 + k*n)/(k*(k-1)/2 + k*m + m*(m-1)/2 + m*n + n*(n-1)/2 + k*n)

def main():
    with open('../data/rosalind_iprb.txt', 'r') as myfile:
        source=myfile.read().replace('\n', '')
    list = [int(x) for x in source.split()]
    print(iprb(list[0], list[1], list[2]))
    #print(iprb(0, 2, 0))
    
if __name__ == '__main__':
    main()