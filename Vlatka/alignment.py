import sys
from time import time

MATCH = 2
MISMATCH = -1
GAP = -2

def match(a, b):
	if a == b:
		return MATCH
	else:
		return MISMATCH

def rev(x):
	return x[::-1]

# returns the position of the cell with maximum value
def SmithWaterman(x, y):
	m = len(x)
	n = len(y)
	h = [[0 for i in range(m+1)], [0]]
	max_value = 0
	max_i = 0
	max_j = 0
	for j in range(1, n+1):
		for i in range(1, m+1):
			h[1].append(max(0, h[0][i-1] + match(x[i-1], y[j-1]), h[1][i-1] + GAP, h[0][i] + GAP))
			if(h[1][i] > max_value):
				max_value = h[1][i]
				max_i = i-1
				max_j = j-1
		h[0] = h[1][:]
		h[1] = [0]
				
	return [max_i, max_j]

# "crops" the initial arrays to smaller sub-arrays that correspond to the optimal local alignment
def LocalToGlobal(x, y):
	border = SmithWaterman(x, y)
	xx = x[:border[0]+1]
	yy = y[:border[1]+1]
	xx_reversed = xx[::-1]
	yy_reversed = yy[::-1]
	border = SmithWaterman(xx_reversed, yy_reversed)
	xxx = xx_reversed[:border[0]+1]
	yyy = yy_reversed[:border[1]+1]
	return [xxx[::-1], yyy[::-1]]
	
# returns the last row of the Needleman-Wunsch matrix	
def NWScore(x, y):
	score = [[0],[]]
	for j in range(1, len(y)+1):
		score[0].append(score[0][j-1] + GAP)
	for i in range(1, len(x)+1):
		score[1].append(score[0][0] + GAP)
		for j in range(1, len(y)+1):
			score[1].append(max(score[0][j-1] + match(x[i-1], y[j-1]), score[1][j-1] + GAP, score[0][j] + GAP))
		score[0] = score[1][:]
		score[1] = []
	return score[0]
	
# finds the optimal global alignmnent computed using Needleman-Wunsch algorithm;
# it's called when one or both arrays has length == 1, so it's using linear space	
def NeedlemanWunsch(x, y):
	f = [[0 for i in range(len(y)+1)] for j in range(len(x)+1)]
	for i in range(len(x)+1):
		f[i][0] = GAP*i
	for j in range(len(y)+1):
		f[0][j] = GAP*j
	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			f[i][j] = max(f[i-1][j-1] + match(x[i-1], y[j-1]), f[i-1][j] + GAP, f[i][j-1]+ GAP)
			
	alx = []
	aly = []
	i = len(x)
	j = len(y)
	while(i > 0 or j > 0):
		if(i > 0 and j > 0 and f[i][j] == f[i-1][j-1] + match(x[i-1], y[j-1])):
			alx.append(x[i-1])
			aly.append(y[j-1])
			i -= 1
			j -= 1
		elif(i > 0 and f[i][j] == f[i-1][j] + GAP):
			alx.append(x[i-1])
			aly.append('-')
			i -= 1
		elif(j > 0 and f[i][j] == f[i][j-1] + GAP):
			alx.append('-')
			aly.append(y[j-1])
			j -= 1
	return ''.join(reversed(aly)), ''.join(reversed(alx))

#finds the index where the sum of elements of two lists is maximal
def PartitionY(sl, sr):
	sum_list = []
	srr = sr[::-1]
	for i in range(len(sl)):
		sum_list.append(sl[i] + srr[i])
	return sum_list.index(max(sum_list))
	
#finds the optimal global alignment of two arrays; linear memory complexity	
def Hirschberg(x, y):
	z = ""
	w = ""
	if len(x) == 0 or len(y) == 0:
		if len(x) == 0:
			z = ""
			w = ""
			for i in range(len(y)):
				z += '-'
				w += y[i]
		elif len(y) == 0:
			for i in range(len(x)):
				z += x[i]
				w += '-'
	elif len(x) == 1 or len(y) == 1:
		z, w = NeedlemanWunsch(x, y)
	else:
		xmid = len(x)/2
		score_l = NWScore(x[:xmid], y)
		score_r = NWScore(rev(x[xmid:]), rev(y))
		ymid = PartitionY(score_l, score_r)
		zl, wl = Hirschberg(x[:xmid], y[:ymid])
		zr, wr = Hirschberg(x[xmid:], y[ymid:])
		z = zl+zr
		w = wl+wr
		
	return z, w
	
def ReadAFile(file_name):
	sequence = open(file_name).read().splitlines()
	if sequence[0][0] == '>':
		sequence.pop(0)
	return ''.join(sequence)
	
def WriteToAFile(file_name, first_seq, sec_seq):
	out_file = open(file_name, 'w')
	out_file.write(">optimal local alignment of the first sequence\n{0}" \
	"\n>optimal local alignment of the second sequence\n{1}".format(first_seq, sec_seq))

if __name__ == "__main__":
	
	if len(sys.argv) < 3:
		sys.exit('Usage: python alignment.py <filename1> <filename2>')
	x = ReadAFile(sys.argv[1])
	y = ReadAFile(sys.argv[2])
	t_start = time()
	cropped_x, cropped_y = LocalToGlobal(x, y)
	aligned_cropped_x, aligned_cropped_y = Hirschberg(cropped_x, cropped_y)
	t_end = time()
	WriteToAFile('output.txt', aligned_cropped_x, aligned_cropped_y)
	print("Time: "+str(t_end-t_start))
	
	
	
	




	
