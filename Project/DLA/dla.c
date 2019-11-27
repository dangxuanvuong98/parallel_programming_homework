#include "mpi.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>

const int numIteration = 1000;
const int log_cycle = 50;
const double epsilon = 0.001;
const double omega = 1.5;
const double eta = 2.0;
const int N = 40;

double max(double a, double b) {
	return a > b ? a : b;
}

double sum(int rank, int size, double *c, int number) {
	double res = 0;
	int i;
	for (i = 0; i < number; i++) {
		res += c[i];
	}
	return res;
}

double _SOR(int rank, int size, double *cNeighbour, double cSelf) {
	double sumNeighbour = sum(cNeighbour, 4);
	return (omega / 4) * (sumNeighbour) + (1 - omega) * cSelf;
}

void updateConcentration(int rank, int size, int rb) {

}

void initialize(int rank, int size, double *concentation, int *bacilli) {
	int i, j;
	for (j = 0; j < N; j++) {
		for (i = 0; i < N - 1; i++) {
			concentation[i][j] = 0;
			bacilli[i][j] = 0;
		}
		c[N - 1][j] = 1;
		bacilli[N - 1][j] = 0;
	}
	bacilli[0][N / 2] = 1;
}

int *devideTask(int rank, int size, double *concentation, int *bacilli) {
	workload = (int *) malloc ((size - 1) * sizeof(int));
	MPI_Request divideTaskRequest;
	int avg_workload = N / (size - 1);
	int over_workload = N % (size - 1);
	workload[0] = 0;
	int r;
	for (r = 1; r <= size - 1; r++) {
		workload[r] = workload[r - 1] + avg_workload;
		if (r <= over_workload) {
			workload[r] += 1;
		}
		int rWorkload - workload[r] - workload[r - 1];
		MPI_Isend(&rWorkload, 1, MPI_INT, r, -1, MPI_COMM_WORLD, &divideTaskRequest);
		MPI_Isend(concentation + N * workload[r - 1], N * rWorkload,
			MPI_DOUBLE, r, -2, MPI_COMM_WORLD, &divideTaskRequest);
		MPI_Isend(bacilli + N * workload[r - 1], N * rWorkload,
			MPI_DOUBLE, r, -3, MPI_COMM_WORLD, &divideTaskRequest)
	}
	return workload;
}

void collectResult(int rank, int size, double *concentation, int *bacilli, int *workload) {
	int i = 0;
	double error, maxError;
	int stop = 0;
	MPI_Status status;
	int r;
	while (true) {
		i++;
		maxError = 0;
		for (r = 1; r <= size; r++) {
			int rWorkload = workload[r] - workload[r - 1];
			MPI_Recv(concentation, N * rWorkload,
				MPI_DOUBLE, r, 1 * i, MPI_COMM_WORLD, &status);
			MPI_Recv(bacilli, N * rWorkload,
				MPI_DOUBLE, r, 2 * i, MPI_COMM_WORLD, &status);
			MPI_Recv(&error, 1, MPI_DOUBLE, r, 3 * i, MPI_COMM_WORLD, &status);
			maxError = max(maxError, error);
		}
		if (maxError < epsilon || i == numIteration) {
			stop = 1;
		}
		for (r = 1; r <= size; r++) {
			MPI_Isend(&stop, 1, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &divideTaskRequest);
		}
		if ((stop == 1) || (i % log_cycle == 0)) {
			char *logFile = (char *) malloc (1024 * sizeof(char));
			sprintf("log_N=%d_eta=%llf_i=%d", N, eta, i);
			FILE *f = fopen(logFile, "w");
			free(logFile);
			int l, m;
			fprintf(f, "%llf\n\n", maxError);
			for (l = 0; l < N; l++) {
				for (m = 0; m < N; m++) {
					fprintf(f, "%d ", concentation[l][m]);
				}
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
			for (l = 0; l < N; l++) {
				for (m = 0; m < N; m++) {
					fprintf(f, "%d ", bacilli[l][m]);
				}
				fprintf(f, "\n");
			}
		}
	}
}

void recvTask(int rank, int size, double *concentation, int *bacilli, int *workload) {
	MPI_Status status;
	MPI_Recv(&workload, 1, MPI_INT, 0, -1, MPI_COMM_WORLD, &status);
	MPI_Recv(concentation, N * workload, MPI_DOUBLE, 0, -2, MPI_COMM_WORLD, &status);
	MPI_Recv(bacilli, N * workload, MPI_DOUBLE, 0, -3, MPI_COMM_WORLD, &status);
}

void main(int argc, char **argv) {

	double *concentation;
	int *bacilli;

	MPI_Init(&argc, &argv);

	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		concentation = (double *) malloc (N * N * sizeof(double));
		bacilli = (int *) malloc (N * N * sizeof(int));
		initialize(rank, size, concentation, bacilli);
		int *workload = devideTask(rank, size, concentation, bacilli);
		collectResult(rank, size, concentation, bacilli, workload);
	} else {

	}

	MPI_Finalize();

	return;
}