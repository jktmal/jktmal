import numpy as np
import uni
import pot
import algorytmy
from itertools import combinations
import pylab as py

def run(steps, dt, syst):
	for k in range(steps):
		algorytmy.leapfrog.update(syst, dt)
	return( np.array([el.X for el in syst.Uni["compounds"]]), [el.V for el in syst.Uni["compounds"]] )

N = 2 # liczba kulek
T_init = 300 
k = 1000 # liczba kroków
dt = 0.01 # długość kroku czasowego

#a = 100
#b = 100
#c = 100

pos = [np.random.random(2)*10-5 for i in range(N)]
vel = [np.random.random(2)*1-0.5 for i in range(N)]
L = [uni.ball(pos[i], vel[i], 0.1) for i in range(N)]
b = {L[i]:[at for at in L[:i]+L[i+1:] ] for i in range(N)}

#print(b)
U = uni.box(L, b)

results = run(k, dt, U)
dist = [np.linalg.norm(results[0][0][i]-results[0][1][i]) for i in range(len(results[0][0]))]
#print(dist)
#print(results[0][0]) #len(results[0]),np.linspace(0, len(results[0][0])))
py.plot(np.linspace(0, len(dist), len(dist)), dist)
py.show()
