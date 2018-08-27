# Generalized Fiboonacci generator with decay
import sys

# n - number of steps, k - multiplication factor 
# dynamic prog with space optimization
def fib(n, k):
    if n == 1 or n == 2:
        return 1
    else:
        fn_1 = 1
        fn_2 = 1
        for i in range (3, n):
            fn = fn_1 + k*fn_2
            fn_2 = fn_1
            fn_1 = fn
        return fn_1 + k*fn_2

# n - number of steps, k - multiplication factor, m - decay factor         
def fibd(n, k , m):
    if n <= m:
        return fib(n, k)
    else:
        f = []
        f.append(1)
        f.append(1)
        for i in range(2, m):
            f.insert(i, f[i-1] + k*f[i-2])
        for i in range(m, n-1):
            fn = k*sum(f[:-1])
            f[:-1] = f[1:]
            f[-1] = fn
        return k*sum(f[:-1])
            
def main():
    with open('../data/rosalind_fibd.txt', 'r') as myfile:
        source=myfile.read().replace('\n', '')
    list = source.split()
    print(fibd(int(list[0]), 1, int(list[1])))
    #print(fibd(92, 1, 17))
    
if __name__ == '__main__':
    main()