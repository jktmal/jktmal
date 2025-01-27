import numpy as np
import uni
import pot
import algorytmy
from itertools import combinations
import pylab as py
import visualisation
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument("-", "--files", nargs='+', type=str, help="nazwy plikow")
parser.add_argument("-n",  type=int, help="number of atoms in crystal's columns") 
parser.add_argument("-k", type=int, help="number of steps in MD simulation")
parser.add_argument("-dx", type=float, help="displacement of chosen column orientated according to X axis ")
parser.add_argument("-j", type=int, help="index of displaced column ")
parser.add_argument("-dt", type=float, help="time step in MD simulation")
parser.add_argument("-T", type=float, help="temp")
#parser.add_argument("-p", type=str, help='path to the folder with results')
#parser.print_help()

args = parser.parse_args()

def run(steps, dt, syst, T):
	#plt.ion()
	#fig = plt.figure()
	#ax = fig.add_subplot(111)
	#data, = ax.plot([el.X[-1][0] for el in syst.Uni["compounds"]], [el.X[-1][1] for el in syst.Uni["compounds"]], '.')
	#plt.ylim((0, syst.d))
	#plt.xlim((0, syst.d))
	#plt.show() tego nie odhaszowywać
	for k in range(steps):
		#visualisation.update_line(fig, data, [el.X[-1] for el in syst.Uni["compounds"]])
		algorytmy.leapfrog.update(syst, dt, T)
		'''for i in range(uni.vneumann.n):
            syst.down[i].X= syst.down[i].X[-1] + syst.Uni['compounds'][n*(n-1)+i].X[-1] - syst.Uni['compounds'][n*(n-1)+i].X[-2] 
            syst.up[i].X= syst.up[i].X[-1] + syst.Uni['compounds'][i].X[-1] - syst.Uni['compounds'][i].X[-2] 
            syst.left[i].X= syst.left[i].X[-1] + syst.Uni['compounds'][n-1+i*n].X[-1] - syst.Uni['compounds'][n-1+i*n].X[-2] 
            syst.right[i].X= syst.right[i].X[-1] + syst.Uni['compounds'][i*n].X[-1] - syst.Uni['compounds'][i*n].X[-2] '''
		print(k)
	#fig.clear()
	#plt.ioff()
	return( np.array([el.X for el in syst.Uni["compounds"]]), [el.V for el in syst.Uni["compounds"]], syst.energy )


N = 2 # liczba kulek
T_init = 300 
k = 1000 # liczba kroków
dt = 0.01 # długość kroku czasowego

#a = 100
#b = 100
#c = 100

#pos = [np.random.random(2)*10-5 for i in range(N)] pierwotne
#vel = [np.random.random(2)*1-0.5 for i in range(N)] pierwotnie

#pos = [np.array([float(i), 0.0]) for i in range(N)] 10 kulek połączonych

#vel = [np.array([0.0, 0.0]) for i in range(N)] 10 kulek połączonych
'''pos[-1] = np.array([pos[-1][0], pos[-1][1]]) # do pierwszego +0.1
L = [uni.ball(pos[i], vel[i], 0.1) for i in range(N)]
#b = {L[i]:[at for at in L[:i]+L[i+1:] ] for i in range(N)}
print([[i, ((int(i))+1)%N] for i in range(N)])
b = {L[i]:[L[(i-1)%N],L[(i+1)%N]]  for i in range(N)}
b[L[0]] = [L[1]]
b[L[-1]] = [L[-2]]'''

#print(pos)
#print(b)
#U = uni.box(L, b)

#U = uni.vneumann(args.n, args.j, args.dx)

U = uni.ljgas(args.n)


results = run(args.k, args.dt, U, args.T)

print(results[2])
np.savetxt('energy.txt', results[2])
plt.plot(np.arange(len(results[2])), results[2])
plt.ylim(np.min(results[2]), np.max(results[2]))
plt.savefig('energy.pdf')
py.show()
#dist = [np.linalg.norm(results[0][0][i]-results[0][-1][i]) for i in range(len(results[0][0]))]
#print(dist)
#print(results[0][0]) #len(results[0]),np.linspace(0, len(results[0][0])))
#py.plot(np.linspace(0, len(dist), len(dist)), dist)
#py.show()
#dist = [results[0][0][i][0] for i in range(len(results[0][0]))]
#dist2 = [results[0][-1][i][0] for i in range(len(results[0][0]))]
#py.plot(np.linspace(0, len(dist), len(dist)), dist, 'b')
#py.plot(np.linspace(0, len(dist), len(dist)), dist2, 'r')
#py.show()

