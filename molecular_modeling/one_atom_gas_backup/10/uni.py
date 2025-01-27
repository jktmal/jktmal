import numpy as np

class box():
	def __init__(self, l, bonds):
		self.Uni = {'compounds': l, 'bonds': bonds}

class vneumann(box):
	n = 0
	def __init__(self, n, k, dx):
		pos = [np.array([float(j), float(i)]) for i in range(n) for j in range(n)]
		vel = [np.array([0.0, 0.0]) for i in range(n*n)]
		L = [ball(pos[i], vel[i], 0.1) for i in range(n*n)]
		#L = [[L[i+j*n] for i in range(n)] for j in range(n)]
		#b = {L[i][j]:[L[i][(j-1)%n], L[(i-1)%n][j], L[i][(j+1)%n], L[(i+1)%n][j]] for i in range(n) for j in range(n)}
		b = {L[i*n+j]:[L[i*n+(j-1)%n], L[((i-1)%n)*n+j], L[i*n+(j+1)%n], L[((i+1)%n)*n+j]] for i in range(n) for j in range(n)}
		
		for i in range(n):
			#print(L[i][k].X)
			L[i*n+k].X = [np.array([L[i*n+k].X[-1][0]+dx, L[i*n+k].X[-1][1]])]
			#L[i*n+k].X = [np.array([L[i*n+k].X[-1][0]+dx, L[i][k].X[-1][1]])]
		down = [ball(pos[i]+np.array([0.0, -1.0]), vel[i], 0.1) for i in range(n)]
		up = [ball(pos[i+n*(n-1)]+np.array([0.0, 1.0]), vel[i+n*(n-1)], 0.1) for i in range(n)]
		left = [ball(pos[i*n]+np.array([-1.0, 0.0]), vel[i*n], 0.1) for i in range(n)]
		right = [ball(pos[n-1+i*n]+np.array([1.0, 0.0]), vel[n-1+i*n], 0.1) for i in range(n)]
		for i in range(1, n-1):
			b[L[i]] = [L[i-1], down[i], L[i+1], L[i+n]]
			b[L[i+n*(n-1)]] = [L[i+n*(n-1)-1], L[i+n*(n-1)-n], L[i+n*(n-1)+1], up[i]]
			b[L[i*n]] = [left[i], L[i*n-n], L[i*n+1], L[i*n+n]]
			b[L[n-1+i*n]] = [L[n-1+i*n-1], L[n-1+i*n-n], right[i], L[n-1+i*n+n]]
		b[L[0]] = [ left[0], down[0], L[1], L[n] ]
		b[L[n-1]] =  [L[n-2], down[-1], right[0], L[2*n-1] ]
		b[L[n*n-n]] = [left[-1], L[n*n-2*n], L[n*n-n+1], up[0] ]
		b[L[n*n-1]] = [L[n*n-2], L[n*n-n-1], right[-1], up[-1] ]
		

		#print([el.X  for l in b.values() for el in l])
		n = n
		self.Uni = {'compounds': L, 'bonds': b}
		self.down = down
		self.up = up
		self.left = left
		self.right = right

class ljgas(box):
	energy = []
	def __init__(self, n):
		pos = [np.array([(n*n)/(n+1)+float(j)*(n*n)/(n+1), (n*n)/(n+1)+float(i)*(n*n)/(n+1)]) for i in range(0, n) for j in range(0, n)] # wcześniej 1 n+1
		vel = [np.array([0.0, 0.0]) for i in range(n*n)]
		L = [ball(pos[i], vel[i], 0.1) for i in range(n*n)]
		n = n
		b = {el:[] for el in L}
		self.Uni = {'compounds': L, 'bonds': b}
		self.d = float(n)*float(n)

class lj2atgas(box):
	def __init__(self, n):
		pos = [np.array([float(j)*(n), float(i)*(n)])+odl for i in range(0, n) for j in range(1, n+1) for odl in [np.array([0.0, -0.5]), np.array([0.0, 0.5])]]
		vel = [np.array([0.0, 0.0]) for i in range(2*n*n)]
		L = [ball(pos[i], vel[i], 0.1) for i in range(2*n*n)]
		n = n
		b = {L[i]:L[i+(-1)**(i%2)] for i in range(len(L))}
		self.Uni = {'compounds': L, 'bonds': b}
		self.d = float(n)


class moore(box):
	n = 0
	def __init__(self, n, k, dx):
		pos = [np.array([float(j), float(i)]) for i in range(n) for j in range(n)]
		vel = [np.array([0.0, 0.0]) for i in range(n*n)]
		L = [ball(pos[i], vel[i], 0.1) for i in range(n*n)]
		#L = [[L[i+j*n] for i in range(n)] for j in range(n)]
		#b = {L[i][j]:[L[i][(j-1)%n], L[(i-1)%n][j], L[i][(j+1)%n], L[(i+1)%n][j]] for i in range(n) for j in range(n)}
		b = {L[i*n+j]:[L[i*n+(j-1)%n], L[((i-1)%n)*n+j], L[i*n+(j+1)%n], L[((i+1)%n)*n+j]] for i in range(n) for j in range(n)}
		
		for i in range(n):
			#print(L[i][k].X)
			L[i*n+k].X = [np.array([L[i*n+k].X[-1][0]+dx, L[i*n+k].X[-1][1]])]
			#L[i*n+k].X = [np.array([L[i*n+k].X[-1][0]+dx, L[i][k].X[-1][1]])]
		down = [ball(pos[i]+np.array([0.0, -1.0]), vel[i], 0.1) for i in range(n)]
		up = [ball(pos[i+n*(n-1)]+np.array([0.0, 1.0]), vel[i+n*(n-1)], 0.1) for i in range(n)]
		left = [ball(pos[i*n]+np.array([-1.0, 0.0]), vel[i*n], 0.1) for i in range(n)]
		right = [ball(pos[n-1+i*n]+np.array([1.0, 0.0]), vel[n-1+i*n], 0.1) for i in range(n)]
		for i in range(1, n-1):
			b[L[i]] = [L[i-1], down[i], L[i+1], L[i+n]]
			b[L[i+n*(n-1)]] = [L[i+n*(n-1)-1], L[i+n*(n-1)-n], L[i+n*(n-1)+1], up[i]]
			b[L[i*n]] = [left[i], L[i*n-n], L[i*n+1], L[i*n+n]]
			b[L[n-1+i*n]] = [L[n-1+i*n-1], L[n-1+i*n-n], right[i], L[n-1+i*n+n]]
		b[L[0]] = [ left[0], down[0], L[1], L[n] ]
		b[L[n-1]] =  [L[n-2], down[-1], right[0], L[2*n-1] ]
		b[L[n*n-n]] = [left[-1], L[n*n-2*n], L[n*n-n+1], up[0] ]
		b[L[n*n-1]] = [L[n*n-2], L[n*n-n-1], right[-1], up[-1] ]
		

		#print([el.X  for l in b.values() for el in l])
		n = n
		self.Uni = {'compounds': L, 'bonds': b}
		self.down = down
		self.up = up
		self.left = left
		self.right = right


class ball():
	def __init__(self, R, V, m):
		self.X = [np.array([R[0], R[1]])]
		self.V = [np.array([V[0], V[1]])]
		self.m = m
		self.f = 0
		self.E = 0
