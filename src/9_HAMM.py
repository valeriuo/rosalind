# Hamming distance between two strings of equal length

def hamm(s, t):
 
    return sum(cs != ct for cs, ct in zip(s,t))
    
def main():
    with open('../data/rosalind_hamm.txt', 'r') as myfile:
        while True:
            s = myfile.readline().replace('\n', '')
            t = myfile.readline().replace('\n', '')
            if not t:
                break            
            if len(t) != len(s):
                continue
        
            print(hamm(s,t))
if __name__ == '__main__':
    main()