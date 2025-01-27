import numpy as np
import argparse
from matplotlib import pylab as py
parser = argparse.ArgumentParser()
parser.add_argument("-n",  type=int, help="number of atoms in crystal's columns")
parser.add_argument("-f",  type=str, help="neighbourhood matrix") 
args = parser.parse_args()

if args.f != None:
	neighbourhood_matrix = np.loadtxt(args.f)

def energia(stany, n, A=args.f):
	global lam
	stany.resize(n, n)
	if A == None:	
		return lam*np.sum([ stany[i,j]*(stany[(i+1)%n, j] + stany[i, (j+1)%n] + stany[(i-1)%n, j] + stany[i, (j-1)%n]) for i in range(n) for j in range(n) ])
	else:
		return lam*np.sum([ stany[i,j]*A[i,k]*stany[i,k] for i in range(n) for j in range(n) for k in range(n) if k!=n] )

A = np.random.random(size=args.n*args.n)
A[A>=0.5] = 1
A[A<1] = -1
A.resize(args.n, args.n)

kroki = 10000 
lam = 1
kt = 0.1
kt = np.linspace(0.01, 10, 100)
energy = np.zeros( (len(kt), kroki+1), dtype='float')

A.resize(1,args.n*args.n)
E_pocz = energia(A, args.n)
energy[0][0] = E_pocz

#kt = [1,0.5,0.2,0.1,0.05,0.02,0.01]

for i in range(len(kt)):
	for j in range(kroki):
		h = np.random.randint(args.n*args.n)
		A.resize(1,args.n*args.n)
		A[0][h] = A[0][h]*(-1)
		E_konc = energia(A, args.n)
		A.resize(1,args.n*args.n)
		if E_konc < E_pocz:
			E_pocz = E_konc
		else:
			if np.random.random()<np.exp(-(E_konc-E_pocz)/kt[i]):
				E_pocz = E_konc
			else:
				A[0][h] = A[0][h]*(-1)
		energy[i, j+1] = E_pocz

sr = np.apply_along_axis( np.mean, 1, energy )
odch = np.apply_along_axis( np.std, 1, energy, ddof=1)

py.plot( kt, sr, 'b.')
py.show()
py.plot( kt, odch, 'b.')
py.show()

#print(E_pocz)	
