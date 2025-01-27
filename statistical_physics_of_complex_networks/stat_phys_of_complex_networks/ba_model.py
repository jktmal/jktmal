# PROJEKT MODELU SIECI EWOLUUJĄCEJ BARABASI-ALBERT 
# wykonał Jan Malinowski na potrzeby przedmiotu Wprowadzenie do fizyki złożoności. Fizyka statystyczna sieci złożonych.
import time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-mode", type=str,help="tryb pracy: sym - symulacja sieci, stat - wczytanie wyników z symulacji oraz analiza statystyczna", required=True)
#parser.add_argument("-N", type=int, help="increase output verbosity")
parser.add_argument("-m_init", type=int, help="liczba wierzchołków w początkowej klice", required=True)
parser.add_argument("-m", type=int, help="liczba krawędzi dołączanych do sieci w każdym kroku", required=True)
parser.add_argument('-t', nargs='+', type=int, help="lista liczb kroków czasowych dla algorytmu w formacie liczba liczba ... np. -t 10 100 1000")
parser.add_argument("-k", type=int,  help="liczba symulacji dla zadanych warunków początkowych")
# OPCJONALNE
parser.add_argument("-o","--output", type=str, help="nazwa prefiksu dla plików do zapisu tekstowego liczności poszczególnych stopni wierzchołków w sieci")
parser.add_argument("-f", "--files", type=str, help="nazwa prefiksu do wczytania plików tekstowych liczności poczególnych stopni wierzchołków w sieci")
parser.add_argument("-w", "--write", type=str, help="nazwa prefiksu dla plików do zapisu wykresu/grafu")
parser.add_argument("-r", "--read", action='store_true', help="wyświetlanie wykresu/grafu")
parser.add_argument("-g", "--graph", action='store_true', help="użycie tej opcji powoduje wyrysowanie grafu, do obejrzenia go konieczne jest użycie również opcji -r")

args = parser.parse_args()
#print(args)
#print(hasattr(args,'o'))
#print(args.o)

if args.mode == 'sym':
	from projekt4b import *
	k = args.k
	ref = 0
	#result = []
	i = k
	#T = [math.floor(el) for el in list(np.linspace(k, k+1500, 10))]
	#Z = np.arange(z+10000000,dtype=int)
	#while l< 3:
	for proba in range(1,k+1):
		for i in args.t:
			a = time.time()
			result, Z, G, hist = symuluj(args.m_init, args.m,i)
			print((time.time()-a, args.m_init, args.m, i))
			if args.output != None:
				#f.write(str(result)+'\n')
				#print(Z)
				np.savetxt(args.output+str(proba)+'_'+str(i)+'.txt', Z, fmt='%1.0f', delimiter=',')
				f = open(args.output+str(proba)+'_'+str(i)+'.txt', 'a')
				f.write(str(result[0])+','+str(result[1])+','+str(result[2])+'\n')
				#ref = result[-1]
				#l += 1
				#i = k
				f.close()
			if args.read == True or args.write != None:
				#print(Z)
				py.plot(Z[0], Z[1]/(np.sum(Z[1]*Z[0])), 'ro')
				py.plot(Z[0], 2*args.m*args.m/(Z[0]*Z[0]*Z[0]), 'b')
				py.title('m_pocz='+str(args.m_init)+', m='+str(args.m)+', t='+str(i)+', i='+str(proba))
				py.xscale('log')
				py.yscale('log')
				#py.errorbar(X, Y, yerr=yer, fmt='o')
				py.xlabel('k')
				py.ylabel('P(k)')
				#py.ylim(-3.5,0)
				if args.write != None:
					py.savefig(args.write+'_'+str(proba)+'_'+str(args.m_init)+'_'+str(args.m)+'_'+str(i)+'p(k)'+'.pdf', format='pdf')
				py.figure()
				py.title('m_pocz='+str(args.m_init)+', m='+str(args.m)+', t='+str(i)+', i='+str(proba))
				#print(hist)	
				for j in range(len(hist)):
					s = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
					py.plot(range(10**j, len(hist[j])+10**j), hist[j], s[j%7])			
				#print(i)
				py.plot(range(1, i), np.array(list(range(1, i)))**(1/2) , 'k--')			
				py.xscale('log')
				py.yscale('log')
				#py.errorbar(X, Z, yerr=zer, fmt='o')
				py.xlabel('t')
				py.ylabel('k(t)')
				#py.ylim(-3.5,0)'''
				if args.write != None:
					py.savefig(args.write+'_'+str(proba)+'_'+str(args.m_init)+'_'+str(args.m)+'_'+str(i)+'k(t)'+'.pdf', format='pdf')
				if args.read == True:
					py.show()
				py.figure()
			if args.graph == True:
				nx.draw_networkx(G, with_labels=False, node_size=100*(100/i), cmap='inferno')
				py.title('m_pocz='+str(args.m_init)+', m='+str(args.m)+', t='+str(i)+', i='+str(proba))
				if args.write != None:
					py.savefig('siec'+args.write+'_'+str(proba)+'_'+str(args.m_init)+'_'+str(args.m)+'_'+str(i)+'.pdf', format='pdf')
				if args.read == True:
					py.show()
				py.figure()
				
#SEKCJA ODNOŚNIE ANALIZY STATYSTYCZNEJ
elif args.mode == 'stat':
	import numpy as np
	from matplotlib import pylab as py
	wynik = [ [] for j in range(len(args.t))]
	srednie = [ [] for j in range(len(args.t))]
	X = np.array(args.t, dtype=int)
	for i in range(1, args.k+1):
		for k in range(len(X)):
			#print(k)
			if 'args.files' != None:
				b = args.files+str(i)+'_'+str(X[k])+'.txt'
			else:
				b = str(i)+'_'+str(X[k])+'.txt'
				
			#f = open(b, 'r')
			U = np.loadtxt(b, delimiter=',', skiprows=2)
			I = np.genfromtxt(b, delimiter=',', skip_footer=1, dtype=int)
			#print(I)
			#print(U)
			wynik[k].append(U[2])
			srednie[k].append(np.sum(I[0]*I[1])/np.sum(I[1]))
			#print(wynik)	
		#print(wynik)
	
	#Y = np.array([np.mean(el) for el in wynik])
	Y = np.array([np.mean(el) for el in wynik])
	yer = np.array([np.std(el, ddof=1) for el in wynik])
	

	#print(Y)
	print(args.write==None)	

	if args.read == True or args.write != None:
		py.plot(X, Y, 'bo')
		py.title('m_pocz='+str(args.m_init)+', m='+str(args.m))
		py.xscale('log')
		py.errorbar(X, Y, yerr=yer, fmt='o')
		py.xlabel('N')
		py.ylabel('alfa')
		py.ylim(-3.5,0)
		py.plot(X, Z, 'ko')
		if args.write != None:
			py.savefig(args.write, format='pdf')
		if args.read == True:
			py.show()
