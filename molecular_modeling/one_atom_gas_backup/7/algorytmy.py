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

	def update(Universe, dt, T): # była 0 = xo
		z = 0
		for at in Universe.Uni['compounds']:
			at.f = np.sum([pot.harmonic.calc_forces(0.1, 1.0, at, part) for part in Universe.Uni['bonds'][at]]) + np.sum([pot.lj.calc_forces(1, 1/1.122, at, part, Universe.d) for part in (Universe.Uni['compounds'][:Universe.Uni['compounds'].index(at)] +  Universe.Uni['compounds'][Universe.Uni['compounds'].index(at)+1:])    ]     ) + pot.langevin.calc_forces(at, 0.75, T)
		#np.sum([el.calc_forces() for el in pot])
			'''print(at.f, z, Universe.Uni['compounds'].index(at))
			input()
			z += 1'''
			at.E = np.sum([pot.harmonic.calc_energy(0.1, 1.0, at, part) for part in Universe.Uni['bonds'][at]]) + np.sum([pot.lj.calc_energy(1, 1/1.122, at, part, Universe.d) for part in (Universe.Uni['compounds'][:Universe.Uni['compounds'].index(at)] +  Universe.Uni['compounds'][Universe.Uni['compounds'].index(at)+1:])    ]     ) + (at.m*np.linalg.norm(at.V[-1])**2)/2
			#print((at.m*np.linalg.norm(at.V[-1])**2)/2)
		Universe.energy.append(np.sum(np.array([at.E for at in Universe.Uni['compounds']])))
		#if Universe.energy[-1] > 10**3:
		#	input()
		for at in Universe.Uni['compounds']:
			at.V.append(at.V[-1] + (at.f/at.m)*dt)
			at.X.append(at.X[-1] + at.V[-1]*dt)
			if (at.X[-1][0]//Universe.d) != 0:
				at.X[-1][0] = at.X[-1][0] + Universe.d*(-1)*(at.X[-1][0]//Universe.d)
			if (at.X[-1][1]//Universe.d) != 0:
				at.X[-1][1] = at.X[-1][1] + Universe.d*(-1)*(at.X[-1][1]//Universe.d)
