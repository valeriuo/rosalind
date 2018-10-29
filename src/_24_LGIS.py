#Longest Increasing Subsequence 
#0 1 2 3 4 5 6 7 8
#8 2 1 6 5 7 4 3 9
import sys

def lms(input, up):
    n = len(input)
    l = []
    ls = []
    d = {}
    v = -1
    for i in range(n):
        m = 0
        for j in range(i):
            if (up and input[i] > input[j] or not up and input[i] < input[j]) and l[j] > m:
                m = l[j]
                v = j
        l.append(m+1)
        d[i] = v

    m = l.index(max(l))
    ls.append(m)
    for i in range(max(l)-1):
        ls.append(d[m])
        m = d[m]
        
    ls.reverse()

    return ls
        

def main():
    with open('../data/rosalind_lgis_test.txt', 'r') as fp:
        n = int(fp.readline())
        perm = [int(x) for x in fp.readline().split()]
        
        print(*[perm[i] for i in lms(perm, True)])
        print(*[perm[i] for i in lms(perm, False)])
    
if __name__ == '__main__':
    main()