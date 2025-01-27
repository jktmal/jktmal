import numpy as np
import argparse
from matplotlib import pylab as py
#import visualisation
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-n",  type=int, help="number of atoms in crystal's columns")
parser.add_argument("-m",  type=int, help="numbera of atoms in crystal's columns")
#parser.add_argument("-f",  type=str, help="neighbourhood matrix") 
args = parser.parse_args()


skos = np.diag_indices(args.n)
skos = ( ( skos[0][:-1], skos[1][1:] ), ( skos[0][1:], skos[1][:-1] ) )
#print(skos)
#if args.f != None:
#	neighbourhood_matrix = np.loadtxt(args.f)

def create_nbmatrix(D, m):
	nb = np.zeros((len(D), len(D)), dtype=int)
	for i in range(len(D)):
		if D[i]-1 in D:
			nb[i, np.where(D==(D[i]-1))[0][0]] = 1
		if D[i]+1 in D:
			nb[i, np.where(D==(D[i]+1))[0][0]] = 1
		if D[i]-m in D:
			nb[i, np.where(D==(D[i]-m))[0][0]] = 1
		if D[i]+m in D:
			nb[i, np.where(D==(D[i]+m))[0][0]] = 1
	return nb

def energia(nbmatrix):
	global skos
	nbmatrix[skos[0][0], skos[0][1]] = 0
	nbmatrix[skos[1][0], skos[1][1]] = 0
	return -np.sum(nbmatrix)

def spin_energy(nbmatrix, spin):
	global skos
	nbmatrix[skos[0][0], skos[0][1]] = 0
	nbmatrix[skos[1][0], skos[1][1]] = 0
	C = np.dot(np.dot(spin, nbmatrix), spin)
	return C


def spin_mc(nbmatrix, spin, kt):
	energy = np.zeros( 3501, dtype='float')
	E_pocz = spin_energy(nbmatrix, spin)
	energy[0] = E_pocz
	for k in range(3500):
		h = np.random.randint(len(spin))
		spin[h] *= -1
		E_konc = spin_energy(nbmatrix, spin)
		if E_konc < E_pocz:
			E_pocz = E_konc
		else:
			if np.random.random()<np.exp(-(E_konc-E_pocz)/kt):
				E_pocz = E_konc
			else:
				spin[h] *= -1
		energy[k+1] = E_pocz
	#print(energy)
	#print(np.mean(energy[len(energy)//2:]))
	return np.mean(energy[len(energy)//2:])

def move_end():
	global A
	if np.random.random()<0.5:
		a = A[1]
		while True:
			if (A[1]-1) in A and (A[1]+args.m) in A and (A[1]+1) in A and (A[1]-args.m) in A:
				return False
			a = np.random.choice([ A[1]-1, A[1]+args.m, A[1]+1, A[1]-args.m ] )
			if a not in A:
				A[0] = a
				break
		return True
	else:
		a = A[-2]
		while True:
			if (A[-2]-1) in A and (A[-2]+args.m) in A and (A[-2]+1) in A and (A[-2]-args.m) in A:
				return False
			a = np.random.choice([ A[-2]-1, A[-2]+args.m, A[-2]+1, A[-2]-args.m ] )
			if a not in A:
				A[-1] = a
				break
		return True
		
def move_square():
	global A
	l = []
	for i in range(1, len(A)-1):
		if (A[i]-A[i-1])!=(A[i+1]-A[i]) and np.abs(A[i+1]-(A[i]-A[i-1])) not in A:
			l.append(i)
	if l != []:
		el = np.random.choice(range(len(l)))
		A[l[el]] = A[l[el]]-(A[l[el]]-A[l[el]-1])+(A[l[el]+1]-A[l[el]])
		return True
	return False

def reptation():
	global A
	#A += A[1]-A[0]
	A[:-1] = A[1:]
	A[-1] = -1
	a = A[-2]
	while True:
		if (A[-2]-1) in A and (A[-2]+args.m) in A and (A[-2]+1) in A and (A[-2]-args.m) in A:
			return False
		a = np.random.choice([ A[-2]-1, A[-2]+args.m, A[-2]+1, A[-2]-args.m ] )
		if a not in A:
			A[-1] = a
			break
	return True

def flip():
	global A
	l = []
	for i in range(1, len(A)-2):
		if np.abs(A[i+1]-A[i]) == np.abs(A[i+2]-A[i-1]) and np.abs(A[i]-A[i-1])!=np.abs(A[i+1]-A[i]) and (A[i]-2*(A[i]-A[i-1])) not in A and (A[i+1]-2*(A[i]-A[i-1])) not in A:
			l.append((i, A[i]-A[i-1]))
	if np.abs(A[1]-A[0])!=np.abs(A[2]-A[1]) and (A[0]+2*(A[2]-A[1])) not in A and (A[1]+2*(A[2]-A[1])) not in A:
		l.append((0, -A[2]+A[1]))
	if np.abs(A[-1]-A[-2])!=np.abs(A[-2]-A[-3]) and (A[-2]-2*(A[-2]-A[-3])) not in A and (A[-1]-2*(A[-2]-A[-3])) not in A:
		l.append((-2, A[-2]-A[-3]))
	
	if l != []:
		el = np.random.choice(range(len(l)))
		A[l[el][0]] = A[l[el][0]]-2*l[el][1]
		A[l[el][0]+1] = A[l[el][0]+1]-2*l[el][1]
		return True
	return False


def spin_gen(n):
	S = np.random.random(size=n)
	S[S<=0.5] = -1
	S[S>0.5] = 1
	return S

movement_methods = [ move_end, move_square, flip, reptation ]
#movement_methods = [ flip ]



S = spin_gen(args.n)

A = np.zeros(args.n, dtype=int)

A[0] = (args.m//2 - args.n//2)*args.m+ args.m//2 - args.n//2
A[1] = np.random.choice([ A[0]-1, A[0]+args.m, A[0]+1, A[0]-args.m ] )

for i in range(2, args.n):
	A[i] = A[i-2]
	if (A[i-1]-1) in A[:i] and (A[i-1]+args.m) in A[:-i] and (A[i-1]+1) in A[:i] and (A[i-1]-args.m) in A[:i]:
		print('tu powinien byc koniec')
		break
	while A[i] in A[:i]:
		A[i] = np.random.choice([ A[i-1]-1, A[i-1]+args.m, A[i-1]+1, A[i-1]-args.m ] )

#A.resize(args.n, args.n)

lock = 0
kroki = 30000 # było 100k
lam = 1
kt = 0.1
kt = np.linspace(1, 8, 20) # wcześniej 0.01, 10, 100
energy = np.zeros( (len(kt), kroki+1), dtype='float')

#A.resize(1,args.n*args.n)
E_pocz = spin_energy(create_nbmatrix(A, args.m), S)
energy[0][0] = E_pocz

#kt = [1,0.5,0.2,0.1,0.05,0.02,0.01]
'''plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
data, = ax.plot(A%args.m, A//args.m, '.')
plt.ylim((0, args.m))
plt.xlim((0, args.m))
'''
for i in range(len(kt)):
	print(kt[i])
	if lock > 100:
		print(str(kt[i-1])+' problem')
	if i%(args.m*args.m)==0 and i!=0 and np.abs(A[args.n//2]-args.m*(args.m//2+1/2)) >= (args.m//4):
		A = A - args.m//4
		print('check')	
	lock = 0
	#input()
	for j in range(kroki):
		#print(energy[i,j])
		print(j)
		accurate_choice = False
		B = np.copy(A)
		while accurate_choice == False:
			A = np.copy(B)
			h = np.random.choice(movement_methods)
			#print(h)
			#A.resize(1,args.n*args.n)
			#A[0][h] = A[0][h]*(-1)
			accurate_choice = h()
		#input()
		S = spin_gen(args.n)
		E_konc = spin_mc(create_nbmatrix(A, args.m), S, kt[i] )
		#A.resize(1,args.n*args.n)
		if E_konc < E_pocz:
			E_pocz = E_konc
			lock = 0
		else:
			if np.random.random()<np.exp(-(E_konc-E_pocz)/kt[i]):
				E_pocz = E_konc
				lock = 0
			else:
				#A[0][h] = A[0][h]*(-1)
				A = np.copy(B)
				lock += 1
		
		energy[i, j+1] = E_pocz
		if lock > 1000:
			break
		#visualisation.update_line(fig, data, (A%args.m, A//args.m))
		#print(A%args.m,A//args.m)

np.savetxt('energia', energy)
'''		
sr = np.apply_along_axis( np.mean, 1, energy )
odch = np.apply_along_axis( np.std, 1, energy, ddof=1)

fig.clear()
plt.ioff()



py.plot( kt, sr, 'b.')
py.show()
py.plot( kt, odch, 'b.')
py.show()
'''
#print(E_pocz)	
