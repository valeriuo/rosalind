#25 - Genome Assembly as Shortest Superstring
# Needleman-Wunch algorithm
import sys

class cell:
	def __init__(self, up, val):
		self.up = up
		self.val = val

def nw(s1, s2, m=1, g=-1, x=-1):
	l1, l2 = len(s1), len(s2)
	score = []
	paths = [[(l2, l1)]]
	best_match = 0
	best_contig = ""
	best_path = []
	o1 = ""
	o2 = ""

	score.append([cell([], 0)] + [cell([(0, i-1)], i*g) for i in range(1,l1+1)])
	for i in range(1,l2+1):
		score.append([cell([], 0) for i in range(l1+1)])
		score[i][0].val = -i
		score[i][0].up = [(i-1, 0)]
		for j in range(1,l1+1):
			d = m if s2[i-1] == s1[j-1] else x
			cand_val = [score[i-1][j-1].val + d, score[i][j-1].val + g, score[i-1][j].val + g]
			cand_ind = [(i-1, j-1), (i,j-1), (i-1, j)]
			score[i][j].val = max(cand_val)
			score[i][j].up += [cand_ind[ind] for ind, mv in enumerate(cand_val) if mv == score[i][j].val]

	# for i in range(l2+1):
		# print(i, end="")
		# for j in range(l1+1):
			# print(*score[i][j].up, end=" | ")
		# print("\n")

	for cpath in paths:
		cpos = cpath[0]
		while cpos != (0, 0):
			l = len(score[cpos[0]][cpos[1]].up)
			if l > 1:
				for i in range(1,l):
					npath = cpath.copy()
					npath.insert(0, score[cpos[0]][cpos[1]].up[i])
					paths.append(npath)
			cpath.insert(0, score[cpos[0]][cpos[1]].up[0])
			cpos = cpath[0]

	for p in paths:
		#print(*p)
		m = 0
		for i in range(1, len(p)):
			if p[i][0] - p[i-1][0] == p[i][1] - p[i-1][1]:
				m += 1
			else:
				if m > best_match:
					best_match = m
					best_path = p
				m = 0

	#print(*best_path)
	for i in range(1, len(best_path)):
		x = best_path[i][0] - best_path[i-1][0]
		y = best_path[i][1] - best_path[i-1][1]
		if x is 1 and y is 1:
			o1 += s1[best_path[i-1][0]]
			o2 += s1[best_path[i-1][0]]
			best_contig += s1[best_path[i-1][1]]
		elif x is 1:
			o1 += "-"
			o2 += s2[best_path[i-1][0]]
			best_contig += s2[best_path[i-1][0]]
		elif y is 1:
			o1 += s1[best_path[i-1][1]]
			o2 += "-"
			best_contig += s1[best_path[i-1][1]]

	#print(o1)
	#print(o2)

	return(best_contig, best_match)

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

	#print(nw("ATTAGACCTGCCGGAA", "GCCGGAATAC", x=-2))
	while len(seq_list) > 1:
		#print(*seq_list)
		best_match = 0
		best_contig = ""
		best_ind = 0
		for i in range(1, len(seq_list)):
			res = nw(seq_list[0], seq_list[i], x=-2)
			if res[1] > best_match:
				best_match = res[1]
				best_contig = res[0]
				best_ind = i
		seq_list.pop(best_ind)
		seq_list.pop(0)
		seq_list.insert(0, best_contig)
	print(seq_list[0])

if __name__ == '__main__':
	main()