import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab as py
import collections
import scipy.optimize as opt
import scipy.linalg
import math
import re
import time
import random

with np.errstate(divide='ignore'):
    np.float64(1.0) / 0.0

#def dodaj_wezel():

def normalizuj(akt,m_pocz,m_dolacz):
        return (1/(m_pocz*(m_pocz-1)+2*m_dolacz*akt), m_pocz*(m_pocz-1)+2*m_dolacz*akt) 

def binary_find(x, L, wsk):
    #print(L[int(len(L)/2)+1:])
    #print(x,L)	
    print(x,L,wsk)
    if x >= L[int(len(L)/2)][0] and x < L[int(len(L)/2)][1]:
        return int(len(L)/2)+wsk
    elif x >= L[int(len(L)/2)][1]:
        return binary_find(x, L[int(len(L)/2)+1:], wsk+int(len(L)/2)+1)	
    elif x < L[int(len(L)/2)][0]:
        return binary_find(x, L[:int(len(L)/2)],wsk) 


def binary_find2(x, L, wsk):
    #print(L[int(len(L)/2)+1:])
    #print(x,L)	
    #print(x,L,wsk)
    #print(x,L,wsk)
    #wait = input("PRESS ENTER TO CONTINUE.")
    if x >= L[int(len(L)/2)-1] and x < L[int(len(L)/2)]:
        return int(len(L)/2)+wsk
    elif x >= L[int(len(L)/2)]:
        return binary_find2(x, L[int(len(L)/2):], wsk+int(len(L)/2))	
    elif x < L[int(len(L)/2)-1]:
        return binary_find2(x, L[:int(len(L)/2)],wsk) 



'''		
m_pocz = 20
m_dolacz = 10
t = 5000
'''
def symuluj(m_pocz, m_dolacz, t):
        G = np.zeros(m_pocz+1+t,dtype=int)
        #print(len(G))
        H = np.array([(m_pocz-1)*k for k in range(m_pocz+1)],dtype=int)
        G[:len(H)] += H        
        #print(G)
        #wait = input("PRESS ENTER TO CONTINUE.")
        m_aktual = m_pocz+1
        for i in range(t):
                #print(G)
                stopien = G[m_aktual-1]
                B = random.sample(range(int(stopien)), m_dolacz)
                G[m_aktual] = G[m_aktual-1]+m_dolacz
                for el in B:
                    G[binary_find2(el,G[:m_aktual],0):m_aktual+1] += 1
                m_aktual += 1
        #print(G)
        '''ROZKLAD'''
        d=collections.Counter()
        for k in range(1,len(G)):
            d[(G[k]-G[k-1])]+=1
        X = np.array(list(d.keys()), dtype=np.float64)
        Y = np.array(list(d.values()), dtype=np.float64)
        #print(X,Y)
        S = np.sum(Y)
        '''py.plot(X, Y/S, 'r.')
        py.yscale('log')
        py.xscale('log')
        py.grid(True)   
        py.show()'''
        wlk = int(len(X)/3+1)
        #print(d.most_common(wlk))
        X2 = np.array([el[0] for el in d.most_common(wlk)])
        Y2 = np.array([el[1] for el in d.most_common(wlk)])

        return((m_dolacz, t, np.polyfit(np.log10(X2), np.log10(Y2/S), 1)[0]),np.array([X,Y], dtype=int))

z = 10
k = 4
l = 4
ref = 0
#result = []

i = k
#T = [math.floor(el) for el in list(np.linspace(k, k+1500, 10))]
#Z = np.arange(z+10000000,dtype=int)

#while l< 3:
for proba in range(1,k+1):
    czasy = [10, 100, 1000, 10000, 100000, 1000000]
    for i in czasy:
        f = open(str(proba)+'_'+str(i), 'w')
        a = time.time()
        result, Z = symuluj(z,l,i)
        print(result,time.time()-a)
        f.write(str(result)+'\n')
        #print(Z)
        np.savetxt(str(proba)+'_'+str(i)+'.txt', Z, fmt='%1.0f', delimiter=',')
        #ref = result[-1]
        #l += 1
        #i = k
        f.close()



#print(d)
#print(result)

'''WYKRES'''
'''
py.plot(X, Y/S, 'r.')
py.yscale('log')
py.xscale('log')
py.grid(True)

#par = np.polyfit(X, Y, 3)
#print(params_final)

#py.plot(X, Y, 'r.')
#py.plot(X, params_final[0]*X**params_final[1], 'b')
#py.plot(X, par[0]*X**3+par[1]*X**2+par[2]*X+par[3], 'b')
py.show()


#nx.draw_networkx(G)'''
