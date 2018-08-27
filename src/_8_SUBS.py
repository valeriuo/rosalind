# Substring indexing

def subs(s, t):
    last = -1
    output = []
    while True:
        last = s.find(t, last+1)
        if last == -1:
            break
        output.append(last+1)
    
    return output
    
def main():
    with open('../data/rosalind_subs.txt', 'r') as myfile:
        while True:
            s = myfile.readline().replace('\n', '')
            t = myfile.readline().replace('\n', '')
            if not t:
                break            
            if len(t) > len(s):
                continue
        
            print(' '.join(map(str, subs(s,t))))
if __name__ == '__main__':
    main()