import numpy as np

class potential():
	def calc_energy():
		pass
	def calc_forces():
		pass
	
class harmonic(potential):
	def calc_energy(k, xo, b1, b2):
		return k/2*(np.linalg.norm(b1.X[-1]-b2.X[-1])-xo)**2
	def calc_forces(k, xo, b1, b2):
		#print(b1)
		#print(b1.X[-1], b2.X[-1], k, xo)
		r=b1.X[-1]-b2.X[-1]
		return -k*(np.linalg.norm(r)-xo)*(r/np.linalg.norm(r))

class lj(potential):
	def calc_energy(eps, ro, b1, b2):
		return 4*eps*( (ro/(np.linalg.norm(b1.X[-1]-b2.X[-1])))**12 - (ro/(np.linalg.norm(b1.X[-1]-b2.X[-1])))**6)
	def calc_forces(eps, ro, b1, b2):
		return 4*eps*( -12*ro**12/(np.linalg.norm(b1.X[-1]-b2.X[-1])**13) + 6*ro**12/(np.linalg.norm(b1.X[-1]-b2.X[-1])**7))		

class langevin(potential):
	def calc_energy():
		pass
	def calc_forces():
		theta = np.random.random()*2*np.pi
		R = np.array([np.cos(theta), np.sin(theta)])*0.1
		return R
