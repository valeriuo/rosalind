# GC proportion
import sys

def gc(input):
    output = 0
    for c in input:
        if c == 'G' or c == 'C':
            output += 1

    return output

def main():
    gc_v = 0
    max_read = ""
    cur_read = ""
    all_v = 0
    max_value = 0
    with open('../data/rosalind_gc.txt', 'r') as fp:
        for line in fp:
            if line[0] == '>':
                cur_read = line.strip(">\n")
                if (all_v != 0):
                    if max_value < 100.0*gc_v/all_v:
                        max_value = 100.0*gc_v/all_v
                        max_read = cur_read
                    gc_v = all_v = 0
            else:
                gc_v += gc(line.strip())
                all_v += len(line.strip())
        if (all_v != 0):
            if max_value < 100.0*gc_v/all_v:
                max_value = 100.0*gc_v/all_v
                max_read = cur_read
    print(max_read) 
    print(max_value)
    
if __name__ == '__main__':
    main()