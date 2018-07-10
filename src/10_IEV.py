# Average number of offsprings with dominant characteristic
    
def main():
    with open('../data/rosalind_iev.txt', 'r') as myfile:
        source=myfile.read().replace('\n', '')
    list = [int(x) for x in source.split()]
    coefs = [1, 1, 1, 0.75, 0.5, 0]
    print(2*sum(x[0] * x[1] for x in zip(list, coefs)))
    
if __name__ == '__main__':
    main()