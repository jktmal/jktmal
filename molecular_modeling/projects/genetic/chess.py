import numpy as np
import itertools as it
import collections

def score(A):
	s = 0
	for i in range(len(A)):
		if (A[(i-1)%len(A)]==0 and A[i]==1 and A[(i+1)%len(A)]==1) or (A[(i-1)%len(A)]==1 and A[i]==0 and A[(i+1)%len(A)]==1):
		#if (A[(i-1)%len(A)]==0 and A[i]==1 and A[(i+1)%len(A)]==1):
			s += 1
	return s

def make_child(A, B):
	sample = np.random.choice(range(8), np.random.randint(8), replace=False)
	C = {}
	x = 0
	for el in it.product(range(2), repeat=3):
		if x in sample:
			C[el] = A[el]
		else:
			C[el] = B[el]
		x += 1

	return C

def generate_rules(l):
	#print(l)
	R = [ (i, rules[i][1]) for i in range(len(l)) ]
	child = make_child(rules[l[0][0]][1], rules[l[1][0]][1])
	mutated_id = np.random.choice([el[0] for el in l])
	tp = rules[mutated_id][1]
	prod = list(it.product(range(2), repeat=3))
	mutated_rule = np.random.randint(len(prod))
	
	if tp[prod[mutated_rule]] == 0:
		tp[prod[mutated_rule]] = 1
	else:
		tp[prod[mutated_rule]] = 0
	extra_rules = [ {el:np.random.randint(2) for el in it.product(range(2), repeat=3)}  for i in range(2) ]
	R.append((4, child))
	R.append((5, tp))
	R.append((6, extra_rules[0]))
	R.append((7, extra_rules[1]))
	return(R)

def evolve(d, B):
	F = np.copy(B)
	#print(d)
	#print(B)
	for j in range(len(B)):
		#print(B)
		#print(B[j])
		#print( d[ (B[(j-1)%len(B)], B[j], B[(j+1)%len(B)])] )
		#print(  (B[(j-1)%len(B)], B[j], B[(j+1)%len(B)] ))
		F[j] = d[ (B[(j-1)%len(B)], B[j], B[(j+1)%len(B)]) ]
	#input()
	return F

N = 20

rules = [ (i, {el:np.random.randint(2) for el in it.product(range(2), repeat=3)})  for i in range(8) ]

S = np.random.random(size=N)
S[S<0.5] = 0
S[S>=0.5] = 1
k = 1000
gen = 10

for j in range(gen):
	scores = collections.Counter()
	best = []
	hist = [[] for i in range(8)]
	#print(rules)
	for el in rules:
		A = np.copy(S)
		hist[el[0]].append(A)
		for i in range(k):
			A = evolve(el[1], A)
			hist[el[0]].append(A)
		scores[el[0]] = score(A)
		best.append(A)
	if j != gen-1:
		rules = generate_rules(scores.most_common(4))
		
print(scores.most_common(1))
print(best[scores.most_common(1)[0][0]])
for el in hist[scores.most_common(1)[0][0]][-10:]:
	print(el)
print(hist[scores.most_common(1)[0][0]][0])

print(rules[scores.most_common(1)[0][0]][1])
#print(hist[el[0]][-10:])
