# Generalized Fiboonacci generator
import sys

# n - number of steps, k - multiplication factor 
def fib(n, k):
    if n == 1:
        return 1
    elif n == 2:
        return 1
    else:
        return k*fib(n-2, k) + fib(n-1, k) 

def main():
    with open(sys.argv[1], 'r') as myfile:
        source=myfile.read().replace('\n', '')
    #source = "5 3"
    list = source.split()
    print(fib(int(list[0]), int(list[1])))
    #print(fib(29, 2))
    
if __name__ == '__main__':
    main()