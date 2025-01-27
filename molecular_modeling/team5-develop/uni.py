import numpy as np

class box():
	def __init__(self, l, bonds):
		self.Uni = {'compounds': l, 'bonds': bonds}

class ball():
	
	def __init__(self, R, V, m):
		self.X = [np.array([R[0], R[1]])]
		self.V = [np.array([V[0], V[1]])]
		self.m = m
		self.f = 0
