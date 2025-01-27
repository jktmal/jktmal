import numpy as np
from matplotlib import pylab as py

wynik = [ [], [], [], [], [], [] ]
X = np.array([10,100,1000,10000,100000,1000000],dtype=int)
for i in range(1,5):
    for k in range(0,6):
        b = str(i)+'_'+str(10**(k+1))
        f = open(b, 'r')
        for line in f:
            wynik[k].append(float(line.lstrip('(').rstrip(')').rstrip('\n').rstrip(')').split(',')[2]))
        f.close()
#print(wynik)

Y = np.array([np.mean(el) for el in wynik])
yer = np.array([np.std(el, ddof=1) for el in wynik])


print(Y)
py.plot(X, Y, 'bo')
py.title('m_pocz=100, m=4')
py.xscale('log')
py.errorbar(X, Y, yerr=yer, fmt='o')
py.xlabel('N')
py.ylabel('alfa')
py.ylim(-3.5,0)
py.show()


#print(A)
#print(A[0,:])#print(A[1,:])
