#include "mpi.h"
#include <stdio.h>

#define M 20
#define Time 0.1
#define dt 0.01
#define dx 0.1
#define D 0.1

void DHB2(float *T, float *dT, int i) {
	int l, r, c;
	c = T[i];
	l = (i == 0) ? 100 : T[i-1];
	r = (i == M-1) ? 25: T[+1];
	dT[i] = D * (l - 2*c + r) / (dx * dx);
}

int main(int argc, char ** argv) {
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int seg = M / size;
	MPI_Status stt;

	float T[M];
	int i;
	
	for (int t = dt; t < Time; t += dt) {
		if (rank == 0) {
			for (i = 0; i < size; i++) {
				MPI_Send(T + i * seg, seg, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
			}
		} else {
			float Ttmp[seg], dTtmp[seg];
			MPI_Recv(Ttmp, seg, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stt);
			for (i = 0; i < seg; i++) {
				DHB2(Ttmp, dTtmp, i);
				Ttmp[i] += dTtmp[i] * dt;
			}
			MPI_Send(Ttmp, seg, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
		}

		if (rank == 0) {
			for (int i = 1; i < size; i++) {
				MPI_Recv(T + i * seg, seg, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &stt);
			}
			for (int i = 0; i < M; i++) {
				printf("%f ", T[i]);
			}
			printf("\n");
		}
	}

	MPI_Finalize();
	return 0;
}
