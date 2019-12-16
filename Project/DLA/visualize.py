import matplotlib as mpl
from matplotlib import pyplot
plt = pyplot
import numpy as np
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

f = open("config.txt", "r")

N = int(f.readline())
# 2 hàm f.readline() tiếp theo để bỏ qua thông tin về trạng thái khởi tạo của trực khuẩn
f.readline()
f.readline()
eta = "{0:.2f}".format(round(float(f.readline()), 2))
omega = "{0:.2f}".format(round(float(f.readline()), 2))
epsilon = "{0:.5f}".format(round(float(f.readline()), 5))
infsm = "{0:.8f}".format(round(float(f.readline()), 8))
numIteration = int(f.readline())
fps = int(f.readline())
skip_frame = int(f.readline())

f.close()

C = np.zeros((N, N))
B = np.zeros((N, N))

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='DLA', artist='Dang Xuan Vuong')
writer = FFMpegWriter(fps=fps, metadata=metadata)

fig = plt.figure()
video_name = "./video/video_N={0:d}_eta={1:s}_omega={2:s}.mp4".format(N, eta, omega)
with writer.saving(fig, video_name, 100):
	for i in range(numIteration):
		if i % skip_frame == 0:
			log_file = "./log/log_N={0:d}_eta={1:s}_omega={2:s}_iter={3:d}.log".format(N, eta, omega, i)
			file = open(log_file, "r")
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
					B[l, m] = float(line[m])
					if B[l, m] == 1:
						C[l, m] = -1

			img = pyplot.imshow(C, interpolation='nearest', cmap='coolwarm', vmin=-1, vmax=1)
			image_name = "./image/log_N={0:d}_eta={1:s}_omega={2:s}_iter={3:d}.jpg".format(N, eta, omega, i)
			pyplot.savefig(image_name)
			print("./image/log_N=" + str(N) + "_eta=" + eta + "_iter=" + str(i) + ".jpg")
			writer.grab_frame()
