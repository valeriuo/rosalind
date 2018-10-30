#Longest Increasing Subsequence 
#0 1 2 3 4 5 6 7 8
#8 2 1 6 5 7 4 3 9
import sys

# Compute longest monotonic sequence
def lms(input, up):
    n = len(input)
    li = [] # array of maximum lengths. l[i] is the length of the sequence ending in input[i]
    d = {}  # d[i] is the index that came before i in the longest sequence containing input[i]

    for i in range(n):
        m = 0
        v = -1
        for j in range(i):
            if (up and input[i] > input[j] or not up and input[i] < input[j]) and li[j] > m:
                m = li[j] # compute the maximum
                v = j     # save the index of the maximum
        li.append(m+1)
        d[i] = v

    m = li.index(max(li))
    ls = [m] # the longest sequence
    for i in range(max(li)-1):
        ls = [d[m]] + ls # go back, starting from the maximum
        m = d[m]

    return ls
        

def main():
    with open('../data/rosalind_lgis.txt', 'r') as fp:
        n = int(fp.readline())
        perm = [int(x) for x in fp.readline().split()]
        
        print(*[perm[i] for i in lms(perm, True)])
        print(*[perm[i] for i in lms(perm, False)])
    
if __name__ == '__main__':
    main()