import numpy as np
import scipy.stats as st

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
		#print(b1, b2)
		#print(r)
		#print(-k*(np.linalg.norm(r)-xo)*(r/np.linalg.norm(r)))
		
		return -k*(np.linalg.norm(r)-xo)*(r/np.linalg.norm(r))

class lj(potential):
	tr = np.array([ [-1.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.0, -1.0], [0.0, 0.0], [-1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1.0] ]) 
	def calc_energy(eps, ro, b1, b2, d):
		cb_pos = lj.tr*d+b2.X[-1]
		#cb_dist = np.array([np.linalg.norm(b1.X[-1]-el) for el in cb_pos])
		cb_dist = np.apply_along_axis(np.linalg.norm, 1, -cb_pos+b1.X[-1])
		min_cb_dist = np.min(cb_dist)
		
		r= lj.tr[np.where(cb_dist==min_cb_dist)[0][0]]*d+b2.X[-1]-b1.X[-1]
		return 4*eps*( (ro/(np.linalg.norm(r)))**12 - (ro/(np.linalg.norm(r)))**6)
	
	def calc_forces(eps, ro, b1, b2, d):
		#print(b1.X, b2.X)
		cb_pos = lj.tr*d+b2.X[-1]
		#print(cb_pos)
		#cb_dist = np.array([np.linalg.norm(b1.X[-1]-el) for el in cb_pos])
		cb_dist = np.apply_along_axis(np.linalg.norm, 1, -cb_pos+b1.X[-1])
		#print(cb_dist)
		min_cb_dist = np.min(cb_dist)
		#print(min_cb_dist)
		#print(lj.tr[np.where(cb_dist==min_cb_dist)][0][0], cb_pos)
		r= lj.tr[np.where(cb_dist==min_cb_dist)[0][0]]*d+b2.X[-1]-b1.X[-1]
		#print(np.linalg.norm(r))
		#print(-4*eps*( -12*ro**12/(np.linalg.norm(r)**13) + 6*ro**12/(np.linalg.norm(r)**7))*(r/np.linalg.norm(r)))
		#input()
		return -4*eps*( -12*ro**12/(np.linalg.norm(r)**13) + 6*ro**12/(np.linalg.norm(r)**7))*(r/np.linalg.norm(r))		

class langevin(potential):
	kb = 860
	def calc_energy():
		pass
	def calc_forces(at, gamma, T):
		#theta = np.random.random()*2*np.pi
		R = -gamma*at.V[-1] + np.sqrt(2*gamma*langevin.kb*T)*st.norm.rvs(size=2)
		#R = np.array([np.cos(theta), np.sin(theta)])*0.1
		return R
