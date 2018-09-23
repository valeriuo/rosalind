#Enumerating Gene Orders (permutations - Heap's algorithm)
import sys

def heap(n, A):
    if n == 1: 
        yield A
    else:
        for i in range(n-1):
            for hp in heap(n-1, A):
                yield hp
            j = 0 if (n % 2) == 1 else i
            A[j],A[n-1] = A[n-1],A[j]
        for hp in heap(n-1, A):
            yield hp

    
def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    
    if n > 0:
        perms = [list(p) for p in heap(n, list(range(n+1))[1:])]
        print(len(perms))
        for p in perms:
            print(*p, sep=' ')
    
if __name__ == '__main__':
    main()