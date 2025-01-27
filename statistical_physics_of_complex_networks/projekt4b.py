import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab as py
import collections
#import scipy.optimize as opt
#import scipy.linalg
#import math
#import re
#import time
import random

with np.errstate(divide='ignore'):
    np.float64(1.0) / 0.0

def binary_find2(x, L, wsk):
	if len(L)==1:
		return L[0]
	#print(L[int(len(L)/2)+1:])
	#print(x,L)
	#print(x,L,wsk)
	#print(x,L,wsk)
	#wait = input("PRESS ENTER TO CONTINUE.")
	# [ 1, 2, 3, 4] 
	# [ 1, 2, 3] OK 
	if x >= L[int(len(L)/2)-1] and x < L[int(len(L)/2)]:
		return int(len(L)/2)-1+wsk #poprawka -1
	elif x >= L[int(len(L)/2)]:
		return binary_find2(x, L[int(len(L)/2):], wsk+int(len(L)/2))	
	elif x < L[int(len(L)/2)-1]:
		return binary_find2(x, L[:int(len(L)/2)], wsk) 


#print(binary_find2(1,[0, 9], 0))

def symuluj(m_pocz, m_dolacz, t):
	G = np.zeros(m_pocz+t,dtype=int) #było +1
	F = nx.complete_graph(m_pocz)
	#print(len(G))
	H = np.array([(m_pocz-1)*k for k in range(m_pocz+1)],dtype=int)
	G[:len(H)] += H        
	#print(G)
	#wait = input("PRESS ENTER TO CONTINUE.")
	m_aktual = m_pocz # było +1
	historia = []
	i = 0
	while 10**i < t:
		historia.append([])
		i += 1
	for i in range(t):
		for k in range(len(historia)):
			#print(i, 10**k)
			if i >= 10**k-1:
				historia[k].append(nx.degree(F, m_pocz+10**k-2))# poprawka dodatkowe -1 
		stopien = G[m_aktual-1]
		#B = random.sample(range(int(stopien)), m_dolacz)
		#print(B)
		F.add_node(m_aktual)
		G[m_aktual] = G[m_aktual-1]+m_dolacz
		znalezione = []
		for j in range(m_dolacz):
			szukany = binary_find2(np.random.randint(int(stopien)),G[:m_aktual],0)
			while szukany in znalezione:
				szukany = binary_find2(np.random.randint(int(stopien)),G[:m_aktual],0)
			znalezione.append(szukany)	
			G[szukany:m_aktual+1] += 1
			F.add_edge(m_aktual,szukany)
			#print(m_aktual,szukany)
		#print(m_aktual)
		m_aktual += 1
	#print(G)
	#print(historia)
	'''ROZKLAD'''
	#historia.reverse()
	d=collections.Counter()
	for k in range(1,len(G)):
		d[(G[k]-G[k-1])]+=1
	X = np.array(list(d.keys()), dtype=np.float64)
	Y = np.array(list(d.values()), dtype=np.float64)
	#print(X,Y)
	S = np.sum(Y)
	wlk = int(len(X)/3+1)
	#print(d.most_common(wlk))
	X2 = np.array([el[0] for el in d.most_common(wlk)])
	Y2 = np.array([el[1] for el in d.most_common(wlk)])
	
	return((m_dolacz, t, np.polyfit(np.log10(X2), np.log10(Y2/S), 1)[0]),np.array([X,Y], dtype=int),F,historia)
