#25 - Genome Assembly as Shortest Superstring
# Needleman-Wunch algorithm
# Smith-Waterman algorithm
# k-mer matching - chosen
import sys
import time
import numpy as np
import re
import cProfile
import itertools as it

def reg_match(s1, s2, k=6):
	l1 = len(s1)
	l2 = len(s2)
	beg = s2[0:k]
	end = s2[l2-k:]
	match = 1
	ind_list = [m.start() for m in re.finditer(beg, s1)]
	for i in range(len(ind_list)):
		start = ind_list[i]
		j=k
		match = 1
		while start+j < l1 and j < l2:
			if s1[start+j] != s2[j]:
				match = 0
				break
			j += 1
		if match == 1:
			contig = s1 + s2[j:]
			return (contig, j, "", "")
	ind_list = [m.start() for m in re.finditer(end, s1)]
	for i in range(len(ind_list)):
		start = ind_list[i]
		j=1
		match = 1
		while start >= j and l2-k >= j:
			if s1[start-j] != s2[l2-k-j]:
				match = 0
				break
			j += 1
		if match == 1:
			contig = s2 + s1[start+k:]
			return (contig, j+k-1, "", "")
	return ("", 0, "", "")

def sw(s1, s2, m=1, x=-1, go=10, gp=2, debug=False):
	l1, l2 = len(s1), len(s2)
	score = np.zeros((l2+1, l1+1), np.int32)
	prev = np.zeros((l2+1, l1+1), np.int32)
	px, py = l1, l2
	best_match = 0
	best_contig = ""
	best_path = []
	o1 = ""
	o2 = ""
	cr = 1
	lf = 2
	up = 4

	prev[0,:] = lf
	prev[:,0] = up
	prev[0,0] = 0

	for i, j in it.product(range(1,l2+1), range(1,l1+1)):
		d = m if s2[i-1] == s1[j-1] else x
		if score[i-1, j-1] + d > score[i, j]:
			score[i, j] = score[i-1, j-1] + d
			prev[i, j] += cr
		if score[i, j-1] - gp > score[i, j]:
			score[i, j] = score[i, j-1] - gp
			prev[i, j] += lf
		if score[i-1, j] - gp > score[i, j]:
			score[i, j] = score[i-1, j] - gp
			prev[i, j] += up
		if score[i, j] > best_match:
			best_match = score[i][j]
			px, py = j, i

	if debug:
		for i in range(l2+1):
			print(i, end=" ")
			for j in range(l1+1):
				print(prev[i][j], end=" | ")
			print("\n")

	if py == l2:
		best_contig = s1[l1:px-1:-1]
		o1 = s1[l1:px-1:-1]
		o2 = "-"*(l1-px)
	elif px == l1:
		best_contig = s2[l2:py-1:-1]
		o1 = "-"*(l2-py)
		o2 = s2[l2:py-1:-1]

	cpos = prev[py, px]
	best_match = 0
	tmp_match = 0
	while cpos > 0:
		pmax = 0
		pdir = 0
		if cpos & 0x1:
			pmax = score[py-1, px-1]
			pdir = cr
		if (cpos>>1) & 0x1:
			if score[py, px-1] > pmax or pdir == 0:
				pmax = score[py, px-1]
				pdir = lf
		if (cpos>>2) & 0x1:
			if score[py, px-1] > pmax or pdir == 0:
				pmax = score[py-1, px]
				pdir = up

		if pdir == cr:
			tmp_match += 1
			if tmp_match > best_match:
				best_match = tmp_match
			best_contig += s1[px-1]
			o1 += s1[px-1]
			o2 += s2[py-1]
			px -= 1
			py -= 1
		elif pdir == up:
			tmp_match = 0
			best_contig += s2[py-1]
			o1 += "-"
			o2 += s2[py-1]
			py -= 1
		elif pdir == lf:
			tmp_match = 0
			best_contig += s1[px-1]
			o1 += s1[px-1]
			o2 += "-"
			px -= 1
		cpos = prev[py, px]

	return(best_contig[::-1], best_match, o1[::-1], o2[::-1])

def nw(s1, s2, m=1, g=-1, x=-1, longest_match=False, debug=False):
	l1, l2 = len(s1), len(s2)
	score = np.zeros((l2+1, l1+1), np.int32)
	prev = np.zeros((l2+1, l1+1), np.int32)
	paths = [[(l2, l1)]]
	best_match = 0
	best_contig = ""
	best_path = []
	o1 = ""
	o2 = ""
	cr = 1
	lf = 2
	up = 4

	prev[0,:] = lf
	prev[:,0] = up
	prev[0,0] = 0

	for j in range(1, l1+1):
		score[0, j] = j*g
	for i in range(1,l2+1):
		score[i, 0] = i*g
		for j in range(1,l1+1):
			d = m if s2[i-1] == s1[j-1] else x
			if score[i-1, j-1] + d > score[i, j]:
				score[i, j] = score[i-1, j-1] + d
				prev[i, j] += cr
			if score[i, j-1] - gp > score[i, j]:
				score[i, j] = score[i, j-1] - gp
				prev[i, j] += lf
			if score[i-1, j] - gp > score[i, j]:
				score[i, j] = score[i-1, j] - gp
				prev[i, j] += up

	if debug:
		for i in range(l2+1):
			print(i, end=" ")
			for j in range(l1+1):
				print(score[i][j], end=" | ")
			print("\n")

	if longest_match:
		for cpath in paths:
			cpos = prev[cpath[0, 0], cpath[0, 1]]
			while cpos != 0:
				if cpos in [3, 7]:
					npath = cpath.copy()
					npath.insert(0, (cpath[0, 0], cpath[0, 1]-1))
					paths.append(npath)
				if cpos in [5, 6, 7]:
					npath = cpath.copy()
					npath.insert(0, (cpath[0, 0]-1, cpath[0, 1]))
					paths.append(npath)
				if cpos in [1, 3, 5, 7]:
					cpath.insert(0, (cpath[0, 0]-1, cpath[0, 1]-1))
				if cpos in [2, 6]:
					cpath.insert(0, (cpath[0, 0], cpath[0, 1]-1))
				if cpos == 4:
					cpath.insert(0, (cpath[0, 0]-1, cpath[0, 1]))
				cpos = prev[cpath[0, 0], cpath[0, 1]]
		if debug:
			print(len(paths))
		for p in paths:
			if debug:
				print(*p)
			ml = 0
			for i in range(1, len(p)):
				if p[i][0] - p[i-1][0] == p[i][1] - p[i-1][1]:
					ml += 1
				else:
					if ml > best_match:
						best_match = ml
						best_path = p
					ml = 0
		if debug:
			print(*best_path)
		for i in range(1, len(best_path)):
			x = best_path[i][0] - best_path[i-1][0]
			y = best_path[i][1] - best_path[i-1][1]
			if x is 1 and y is 1:
				o1 += s1[best_path[i-1][1]]
				o2 += s2[best_path[i-1][0]]
				best_contig += s1[best_path[i-1][1]]
			elif x is 1:
				o1 += "-"
				o2 += s2[best_path[i-1][0]]
				best_contig += s2[best_path[i-1][0]]
			elif y is 1:
				o1 += s1[best_path[i-1][1]]
				o2 += "-"
				best_contig += s1[best_path[i-1][1]]
	else: #longest_match
		x, y = l1, l2
		cpos = prev[y, x]
		tmp_match = 0
		while cpos > 0:
			if cpos in [1, 3, 5, 7]:
				tmp_match += 1
				if tmp_match > best_match:
					best_match = tmp_match
				best_contig += s1[x-1]
				o1 += s1[x-1]
				o2 += s2[y-1]
				x -= 1
				y -= 1
			if cpos in [4, 6]:
				tmp_match = 0
				best_contig += s2[y-1]
				o1 += "-"
				o2 += s2[y-1]
				y -= 1
			if cpos == 2:
				tmp_match = 0
				best_contig += s1[x-1]
				o1 += s1[x-1]
				o2 += "-"
				x -= 1
			cpos = prev[y, x]

	return(best_contig[::-1], best_match, o1[::-1], o2[::-1])

def main():
	seq_list = []
	#with open('../data/test_long2.txt', 'r') as f:
	with open('../data/rosalind_long.txt', 'r') as f:
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

	#bs1 = "CTTCCCTAACGTAGAAGTCGCTCACCGTTGGTTGTTCTCTTCGCACGCTACCCGCACACTATATCAGCCAGCACCACTCGCTTCCTATGCACGGATATTATATGCAGCGGCACTTAGCGTACAGCAATCAGAGTCCCAAATCTTTAGGACATAACACGAATTCGGTGACACTACCCACGCGGATAATGCAGCGTGGGATAATCTTTCAGCCACACGCCCTAACTAAACAAAAGTCTCCCTCGATCGCTCAGCGTGGGTGAACATAGTCAATAGCAAATAAGGGCCCGCAGCAGGTCTGCAACGGTTTGGAGCTGTCCGCTCAACCCCAAATCTGCACCTCAGATGGCTCCTTCGAGCCGTGCACTAGAGTACCAATTCCTATCATTACCGTCACATAATACTCCGGCCCGAATTCTGAAATGCTTCTGCGGGGTTTACTTACTGTGTTTAGAGAGGTATCAAGCCTCAATGCATGATGAGAATGTCACACTTCCCTTATGACTCCAAATTCTATAGTACCCGTTTTCACGTTGCACTATGGCTTTTCTATTTATACCTTAAGTTACTGTTAAAAACATACGGTCCCAAATCAGGGAGCTCAATCCCTCGTCGTGTTCGGGTTGCGCATCTGTCTAGGAGGGGTAGCAATAACTCTGCTGTGTGGTACAAGAAGGAAAAGACAGGGTGTTTGCAAACCCCCCGAGTGAGACGTGCACCCCACAAATAGGTGGTGGGAATGTATTCCCCGTGGAGCCAGTCCATCGGTGAGGTCCTCAAATGCAGTGAAACTATGTCTTGAGGAAGACGGGCGGAATCCCTAACCTCGCTAGGCATGCACCGTAGTCATAGATGCGAGTAATTTGATCATCCTCAGCGAAACACAGCGACTCCAATGTAGCTATACATGATGTGTGTCACTCCTTTTGCATACGACGGTCGCTCACAGTAACGTTAGGTGGATAACCCACTGATATGACGATATTACCGTATGGTATTGACCGAGAGAGCCGGAAACCTCAAGTATTCCTAGCGCGATGTGAGCTAATATACACACTGTCGTACTATATGCTTGCCCTTACAGACGCTTGAGGCCGCTACAGAGTGGGACTGTGTCGTAGTGTTACACGAGTATCACTCTCTGAGCGGTTATCATTTTTACTTCAATACTATCAATATACGCCAATTGGGTTGTACAGTGCGAGTCAGAACATACTCGCGCAAACGCGAGCGCCGCGCACTGCCAGCATGTTGTGCTTCCATAGTGCTCGGGTTATACTGTCGTTTAATATCATTGGTTAGCGCAAAGACGGGCATATGATGGCGAAATTCTATAAGTGCCACAACCAAAAGGACGGATCACCAAGTTTAAGAGGCTCGCCTATAGTGAGTCTATCTCCGAATCTTTTAATCAAGTGCGGGCTCGTGTTGTCCCGGCGTAACGCGTGTGAGGAATCCATTACATTATTTGTTAGTCCAGAATAGCCGTTAGAAAGCTAGGGCAATAGGGGCGTCGCTCCAGAGTGCGAATGAAGTTTCACTACACTACCTTCGGTAGGGCCCTGCGTTCTCGAACAGTCCCACTGCGGCGTCCGGTTTTAGGCGTATTAGAGATTATAGATACCGACTCGAAGTAAATAGCTCCGGTGAACTACTGCTTGATAACGTGCTATGGGGTCCGGATACATTATATTAGCGTAGACTCGCCGGCATAAGATCTCATTTACAATTGACCTGAAGCGGCCGCCACGTAACAGCCACTAAGGCCGACTTAATGGCTTGGCACGCCGTCCCTAGCCGGTGTTAAGGTACCGCCCACCCATCTCTTTAATTTTAGTTCCTTAGAGCCCTCGGGGCCCCTAAAGGGTCAGCAATCCTTTTCACCTAATCGTTATTGACTCACATACTTCCCTGTGGTATTGGTACTCCAGCGTTCGGTGCAATATGACCCGATGGCTCATGTGCCGGGTCACCCGTGCCCGACAGCAACTGAGGAGAGAGTAACAGTCACTATAAAGTGTGACAAGTTTACTAGTCCGTGCGGGTTGTCCAAACCGCTGTGGGTGGCTATAGGTTCGGTTACGACAAGGCTTGGGTCTTCTCGTGCTGCTATTCTTATCCTCTGTTCATTGTGTCGAACGTACAAGGGTGATACAGGCGCGTCCTTCATCTTTGACCATATTCAGCGCTAAACTGCTGGCACTGCAGACTGTGAACTCAAGCAAAGTTGTGACTCGAGCGCGCCTCACCTCCTAGGGCTGTCTATCTGTCCATTTCACGTTGGCGTCCCGCCAAGGTTTCTCAGTAGATATATCAGCCTATCAACCATTGTCTTCCTGCTACAGAACGTGCGCGGGCCTTAAGTACCGGGCTCCACTACAGTCATCCACTTTTCAGAGCCGGTCATTCGGAGGAGCGCCAGCTCCGCTATACCGGGGCGAGGGCTGAGTACCCCGAACTCGTAAGATGTGCACATCCATTTTAAACGGTATTCTATGCGGATTAGTCTACATACGTGGCCCCTTGGGTTGGGGTGTATTATCTAATCAAGCGCGTGCCGATGTACAATTACTCCGTCCCCTCCGGGCGTATTAAACGTTCAACCAGCTACGTATCGATCCGGCTTCCGTTTCCGTGCTGGGGTAAGCCTTGTCTCGGGCAATGTTTTGAGACACACGATAAACGCAAGAAGCGAGTCACGCAGTTCAGAGACGGTCTTACGTGAACCTTTCTGTGTACTTGCCGGAGTGGGTCCAAAATGCGCTGCGTGTTTCAATAGTGGCCCCAAGGGGTTGGGAACCCAAGATAGACTTCCCCCGCCATGTTGACCCCGCGTGAAAGCACTTGGAGGGCTTCGGACATGACGAGGGTTCTAGTAGCGGCTCTGAACTAGTAAGCGCATTATCCTTAGCATATAACGCGTATCTTGCGATTACTATTCCCAAACCGCCAAATGTTGGGAGAGGAGGGTAGAGGTCGTGTTCACCTTGCGGGAGAGCATATAACTATAAACCGGCGTAGGTTCGACTACTTCCAGACTGTGCTCTTGATCAACGATCACCAACCCGCGCCTCACTGGTTTCTTTAGTCATCGCCCCGATATAATAATGATGGTCCCCGCCAAGGATGTGTCATACAAAGGTGATTCCTCAACCAGGTAATTGTCACGGACCAAGACTGCACCAGCTGGCGAACAGGTTGGCTATAAATGACCGATTGGTATGCCCCAGCTTTAGTAGGGTAACGTAGCTGGATGCCGCTCGCCGCGGATGAAAGGACGGGATACCCTCACGGAGTTAGGTTGTAGGGTCCTGAGAGTTCATCGGACGTCTCATATTACATGACCATCGCCCCTTCATAGCACGATATTAATTCATGCACTTGGCGATTATGATAGTGTATCATATACTAGCAGGGGTGCGGACATGTGAGGCGATTAAAGTAGTTCATCCAGAAAACAAAGCTAGCGGTGAAATAACTTCACCACTCTTGGCCGAGCTTAGAAATCCCGATGAGGGTACGATGCAGCCAGCCTCGGTCACGCTACAAGCGGGTATTATTCGTTTAAATGCCTATTCTAATGATATTGAAAGGAGTCAGCGTATAACATCGGAGTGACCATATGGAAACTCTCTGATGCTGCTTCTAAGAAGTCGCCGAATCATGTTTCTATTCAAAGGAACTTTATATAATGTGGTCAAGTACCGTAAAGGTTAGTCTCGTCGTATGGGATTCGCGACGGTCAGACATTAGTGTTTTAATAGAAACGTATCAGATTAAACTTTACGCCGTTGTGA"
	#bs2 = "AGACTGCACCAGCTGGCGAACAGGTTGGCTATAAATGACCGATTGGTATGCCCCAGCTTTAGTAGGGTAACGTAGCTGGATGCCGCTCGCCGCGGATGAAAGGACGGGATACCCTCACGGAGTTAGGTTGTAGGGTCCTGAGAGTTCATCGGACGTCTCATATTACATGACCATCGCCCCTTCATAGCACGATATTAATTCATGCACTTGGCGATTATGATAGTGTATCATATACTAGCAGGGGTGCGGACATGTGAGGCGATTAAAGTAGTTCATCCAGAAAACAAAGCTAGCGGTGAAATAACTTCACCACTCTTGGCCGAGCTTAGAAATCCCGATGAGGGTACGATGCAGCCAGCCTCGGTCACGCTACAAGCGGGTATTATTCGTTTAAATGCCTATTCTAATGATATTGAAAGGAGTCAGCGTATAACATCGGAGTGACCATATGGAAACTCTCTGATGCTGCTTCTAAGAAGTCGCCGAATCATGTTTCTATTCAAAGGAACTTTATATAATGTGGTCAAGTACCGTAAAGGTTAGTCTCGTCGTATGGGATTCGCGACGGTCAGACATTAGTGTTTTAATAGAAACGTATCAGATTAAACTTTACGCCGTTGTGATCCTGGATGACAAGTGCTTAACAAAGGGTTTTGTTTGGGCCCAACACCTCGCGTTGGGTTTCGGCTCCATACCTGATCTTACGGGACGGTGGGCGCCCGGCTATCGTATGTAAGCATCGAAGCTACTACGAACTCATGAGCCCAGGCTACTGGCCTAAGTGCCAGTATACGATCTGCTTTGTGCGCTGAACGTGTGGCTATAGTAGGCACTGTCCCTGAGAATTCGTTGATTCTATGGCCGTGCCTGACAACCATACATTTCCTTTGTAAATTCAACCAAGGACCATGCACACTTGTAGTGCAATATTATTGTGGTTAAAGTAAACGAACCACTGCGTACTAGAAGGGCGTCGTGGTGACCCAGGCG"

	#bs1 = "CATATGGTGACCTTAGCCTTGTGATAGCACCTTAGGTCGCTAGACACCTTAG"
	#bs2 = "GACTAGCATATGGTGACCTTAG"

	#t0 = time.process_time()
	#res = sw(bs1, bs2)
	#t1 = time.process_time()
	#print("contig=", res[0], "len=", res[1], "sw done in:", t1 - t0)
	#res = reg_match(bs1, bs2)
	#t2 = time.process_time()
	#print("contig=", res[0], "len=", res[1], "rm done in:", t2 - t1)
	
	#start, end = smith_waterman(bs1, bs2)
	#print(start, end)
	#print(bs1[start:end])
	#print(nw("ATTAGACCTG", "AGACCTGCCG", debug=False))
	#print(sw("ATTAGACCTG", "AGACCTGCCG", debug=False))

	seq_threshold = len(seq_list[0])/10
	while len(seq_list) > 1:
		best_contig = ""
		best_ind = 0
		for i in range(1, len(seq_list)):
			#res = sw(seq_list[0], seq_list[i])
			res = reg_match(seq_list[0], seq_list[i])
			if res[1] >= seq_threshold:
				best_contig = res[0]
				#print(i, res[1])
				#print(best_contig)
				best_ind = i
				break
		seq_list.pop(best_ind)
		seq_list.pop(0)
		seq_list.insert(0, best_contig)
	print(seq_list[0])

if __name__ == '__main__':
	#cProfile.run("main()")
	main()
