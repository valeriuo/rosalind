# Mendel's Second Law - Independent Alleles
# http://rosalind.info/problems/lia/
import sys
import operator as op
import math

def combinations(n, k):
    # k = min(k, n-k)
    # numer = reduce(op.mul, xrange(n, n-k, -1), 1)
    # denom = reduce(op.mul, xrange(1, k+1), 1)
    # return numer//denom
    
    return math.factorial(n)/math.factorial(k)/math.factorial(n-k)
    

def main():
    with open('../data/rosalind_lia.txt', 'r') as myfile:
        source=myfile.read().replace('\n', '')
    list = [int(x) for x in source.split()]
    r = 2**list[0] # number of offsprings on k level
    N = list[1]
    
    print("%0.3f" % (sum([3**(r-i)*combinations(r,i) for i in range(N, r+1)])/4**r))
    
if __name__ == '__main__':
    main()