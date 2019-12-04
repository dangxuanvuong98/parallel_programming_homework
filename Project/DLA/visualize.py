import matplotlib as mpl
from matplotlib import pyplot
import numpy as np

numIteration = 1000
epsilon = 1e-3
eps = 1e-10
omega = 1.5
eta = 1.0
N = 40

C = np.zeros((N, N))
B = np.zeros((N, N))

for i in range(numIteration):
	if i % 100 == 0:
		file_name = "./log/log_N=" + str(N) + "_eta=1.000000_iter=" + str(i) + ".log"
		file = open(file_name, "r")
		for l in range(N):
			line = file.readline()
			line = line.split(' ')
			for m in range(N):
				C[l, m] = float(line[m])

		file.readline()

		for l in range(N):
			line = file.readline()
			line = line.split(' ')
			for m in range(N):
				B[l, m] = int(line[m])

		img1 = pyplot.imshow(C, interpolation='nearest', cmap='coolwarm')
		img2 = pyplot.imshow(B, cmap='Greys')
		pyplot.show()