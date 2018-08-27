# Generalized Fiboonacci generator
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

# recursive version        
def fib_(n, k):
    if n == 1 or n == 2:
        return 1
    else:
        return fib(n-1, k) + k*fib(n-2, k)

def main():
    with open('../data/rosalind_fib.txt', 'r') as myfile:
        source=myfile.read().replace('\n', '')
    #source = "5 3"
    list = source.split()
    print(fib(int(list[0]), int(list[1])))
    #print(fib_(53, 5))
    #print(fib(53, 5))
    
if __name__ == '__main__':
    main()