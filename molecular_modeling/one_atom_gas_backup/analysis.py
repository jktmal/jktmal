import numpy as np
from matplotlib import pylab as py

A = np.zeros((5000, 10), dtype='float')
for i in range(1, 11):
	A[:, i-1] = np.loadtxt(str(i)+'/energy.txt')

sr = np.zeros(10, dtype='float')
odch = np.zeros(10, dtype='float')
T = np.linspace(0.01, 10, 10)

sr[0] = np.mean(A[2500:, 0])
odch[0] = np.std(A[2500:, 0], ddof=1)

sr[1] = np.mean(A[2000:, 1])
odch[1] = np.std(A[2000:, 1], ddof=1)

sr[2] = np.mean(A[3800:, 2])
odch[2] = np.std(A[3800:, 2], ddof=1)
sr[3] = np.mean(A[1700:, 3])
odch[3] = np.std(A[1700:, 3], ddof=1)
sr[4] = np.mean(A[1500:, 4])
odch[4] = np.std(A[1500:, 4], ddof=1)
sr[5] = np.mean(A[1300:, 5])
odch[5] = np.std(A[1300:, 5], ddof=1)
sr[6] = np.mean(A[2300:, 6])
odch[6] = np.std(A[2300:, 6], ddof=1)
sr[7] = np.mean(A[1500:, 7])
odch[7] = np.std(A[1500:, 7], ddof=1)
sr[8] = np.mean(A[3800:, 8])
odch[8] = np.std(A[3800:, 8], ddof=1)
sr[9] = np.mean(A[1300:, 9])
odch[9] = np.std(A[1300:, 9], ddof=1)

py.plot(T, sr, 'bo')
py.show()
py.plot(T, odch, 'bo')
py.show()
