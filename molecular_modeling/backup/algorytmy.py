import numpy as np
import pot

class algorithm():
	def update():
		pass

class leapfrog(algorithm):
	#def __init__(self, dt, syst, pot):
	#	self.dt = dt
	#	self.syst = syst
	#	self.pot = pot

	def update(Universe, dt):
		for at in Universe.Uni['compounds']:
			at.f = np.sum([pot.harmonic.calc_forces(0.1, 0, at, part) for part in Universe.Uni['bonds'][at]]) +pot.langevin.calc_forces()
			       #np.sum([langevin.calc_forces(at, part) for part in Universe['bonds'][at]]) 
		#np.sum([el.calc_forces() for el in pot])	
		
		for at in Universe.Uni['compounds']:
			at.V.append(at.V[-1] + (at.f/at.m)*dt)
			at.X.append(at.X[-1] + at.V[-1]*dt)
			

