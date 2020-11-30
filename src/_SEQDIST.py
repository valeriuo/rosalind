#25 - Genome Assembly as Shortest Superstring
import sys
import numpy as np

def seq_dist(s1, s2):
	l1, l2 = len(s1), len(s2)
	match = []
	dist = 0

	match.append([i for i in range(l1+1)])
	for i in range(1,l2+1):
		match.append([0 for i in range(l1+1)])
		match[i][0] = i
		for j in range(1,l1+1):
			match[i][j] = min(match[i-1][j-1], match[i-1][j], match[i][j-1]) + (0 if s2[i-1] == s1[j-1] else 1)
			
	d = 0
	s = ''
	i, j = l2, l1
	while i > 0 and j > 0:
		if i is 0 or match[i-1][j-1] > match[i][j-1]:
			if d > dist:
				dist = d
			d = 0
			s = s1[j-1]+s
			print("INS: {}, {}, {}, {}".format(i,j,s1[j-1],s))
			j = j-1
		elif j is 0 or match[i-1][j-1] > match[i-1][j] or j == 0:
			if d > dist:
				dist = d
			d = 0
			s = s2[i-1]+s
			print("DEL: {}, {}, {}, {}".format(i,j,s2[i-1],s))
			i = i-1
		else:
			d = d+1
			s = s2[i-1]+s
			print("MAT: {}, {}, {}, {}".format(i,j,s2[i-1],s))
			i = i-1
			j = j-1

	return (s, dist)

def main():
	seq_list = []
	with open('../data/test_long.txt', 'r') as f:
		seq = ''
		for line in f:
			if line[0] != '>':
				seq += line.replace('\n', '')
			else:
				if len(seq) > 0:
					seq_list.append(seq)
					seq = ''
		if	len(seq) > 0:
			seq_list.append(seq)

	print(seq_dist(seq_list[0], seq_list[2]))
	
if __name__ == '__main__':
	main()