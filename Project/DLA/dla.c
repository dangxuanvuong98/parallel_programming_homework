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
const int ran = 1e6;

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

double _SOR(double *cNeighbour, double cSelf) {
	double sumNeighbour = sum(cNeighbour, 4);
	return (omega / 4) * (sumNeighbour) + (1 - omega) * cSelf;
}

void initialize(int rank, int size, double *concentration, int *bacilli) {
	int l, m;
	for (m = 0; m < N; m++) {
		for (l = 0; l < N - 1; l++) {
			concentration[N * l + m] = 0;
			bacilli[N * l + m] = 0;
		}
		c[N * (N - 1) + j] = 1;
		bacilli[N * (N - 1) + j] = 0;
	}
	bacilli[N / 2] = 1;
}

int *devideTask(int rank, int size, double *concentration, int *bacilli) {
	workload = (int *) malloc (size * sizeof(int));
	MPI_Request divideTaskRequest;
	int avg_workload = N / (size - 1);
	int over_workload = N % (size - 1);
	workload[0] = 0;
	int r;
	int stop = 0;
	for (r = 1; r <= size - 1; r++) {
		workload[r] = workload[r - 1] + avg_workload;
		if (r <= over_workload) {
			workload[r] += 1;
		}
		int rWorkload = workload[r] - workload[r - 1];
		MPI_Isend(workload, size, MPI_INT, r, -1, MPI_COMM_WORLD, &divideTaskRequest);
		MPI_Isend(concentration + N * workload[r - 1], N * rWorkload,
			MPI_DOUBLE, r, -2, MPI_COMM_WORLD, &divideTaskRequest);
		MPI_Isend(bacilli + N * workload[r - 1], N * rWorkload,
			MPI_INT, r, -3, MPI_COMM_WORLD, &divideTaskRequest);
		MPI_Isend(&stop, 1, MPI_INT, r, 0, MPI_COMM_WORLD, &divideTaskRequest);
	}
	return workload;
}

void collectResult(int rank, int size, double *concentration, int *bacilli, int *workload) {
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
			MPI_Recv(concentration, N * rWorkload,
				MPI_DOUBLE, r, 1 * i, MPI_COMM_WORLD, &status);
			MPI_Recv(bacilli, N * rWorkload,
				MPI_INT, r, 2 * i, MPI_COMM_WORLD, &status);
			MPI_Recv(&error, 1, MPI_DOUBLE, r, 3 * i, MPI_COMM_WORLD, &status);
			maxError = max(maxError, error);
		}
		if (maxError < epsilon) {
			stop = 1;
		}
		for (r = 1; r < size; r++) {
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
					fprintf(f, "%d ", concentration[N * l + m]);
				}
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
			for (l = 0; l < N; l++) {
				for (m = 0; m < N; m++) {
					fprintf(f, "%d ", bacilli[N * l + m]);
				}
				fprintf(f, "\n");
			}
		}
	}
}

void recvTask(int rank, int size, double *concentration, int *bacilli, int *workload) {
	MPI_Status status;
	workload = (int *) malloc (size * sizeof(int));
	MPI_Recv(workload, size, MPI_INT, 0, -1, MPI_COMM_WORLD, &status);
	concentration = (double *) malloc (N * workload * sizeof(double));
	MPI_Recv(concentration, N * workload, MPI_DOUBLE, 0, -2, MPI_COMM_WORLD, &status);
	bacilli = (int *) malloc (N * workload * sizeof(int));
	MPI_Recv(bacilli, N * workload, MPI_INT, 0, -3, MPI_COMM_WORLD, &status);
}

int isCandidate() {
	return 1;
}

double updateConcentration(int rank, int size, int rb, int *workload,
	double *concentration, int *bacilli
	double *aboveConcentration, double *belowConcentration,
	int *aboveBacilli, double *aboveBacilli) {
	double maxError = 0;
	if (rank > 1) {
		MPI_Isend(concentration, N, MPI_DOUBLE, rank - 1, i,
			MPI_COMM_WORLD, &sendBoundaryRequest);
		MPI_Isend(bacilli, N, MPI_INT, rank - 1, i,
			MPI_COMM_WORLD, &sendBoundaryRequest);
	}
	if (rank < size - 1) {
		MPI_Isend(concentration + N * (workload - 1), N, MPI_DOUBLE, rank + 1, i,
			MPI_COMM_WORLD, &sendBoundaryRequest);
		MPI_Isend(bacilli + N * (workload - 1), N, MPI_INT, rank + 1, i,
			MPI_COMM_WORLD, &sendBoundaryRequest);
	}
	if (rank > 1) {
		MPI_Recv(belowConcentration, N, MPI_DOUBLE, rank - 1, i,
			MPI_COMM_WORLD, &status);
		MPI_Recv(belowBacilli, N, MPI_INT, rank - 1, i,
			MPI_COMM_WORLD, &status);
	}
	if (rank < size - 1) {
		MPI_Recv(aboveConcentration, N, MPI_DOUBLE, rank + 1, i,
			MPI_COMM_WORLD, &status);
		MPI_Recv(aboveBacilli, N, MPI_INT, rank + 1, i,
			MPI_COMM_WORLD, &status);
	}

	int l, m;
	int rWorkload = workload[rank] - workload[rank - 1];
	for (l = 0; l < rWorkload; l++) {
		for (m = 0; m < N; m++) {
			int tmp = l + r + N * workload[rank - 1];
			if (tmp % 2 == rb) {
				double cNeighbour[4];
				if (l == 0) {
					cNeighbour[0] = aboveConcentration[m];
				} else {
					cNeighbour[0] = concentration[l - 1][m];
				}
				if (l == N - 1) {
					cNeighbour[1] = belowConcentration[m];
				} else {
					cNeighbour[1] = concentration[l + 1][m];
				}
				if (m == 0) {
					cNeighbour[2] = concentration[l][N - 1];
				} else {
					cNeighbour[2] = concentration[l][m - 1];
				}
				if (m == N - 1) {
					cNeighbour[3] = concentration[l][0];
				} else {
					cNeighbour[3] = concentration[l][m + 1];
				}
				double oldConcentration = concentration[l][m];
				concentration[l][m] = _SOR(cNeighbour, concentration[l][m]);
				maxError = max(maxError, abs(concentration[l][m] - oldConcentration));
			}
		}
	}
	return maxError;
}

double sumConcentration(double *concentration, int rWorkload) {
	int l, m;
	double res = 0;
	for (l = 0; l < rWorkload; l++) {
		for (m = 0; m < N; m++) {
			res += pow(concentration[l][m], eta);
		}
	}
	return res;
}

void solveTask(int rank, int size, double *concentration, int *bacilli, int *workload) {
	int i = 0;
	int stop;
	MPI_Status status;
	MPI_Request sendBoundaryRequest;
	MPI_Request reportTaskRequest;
	double *aboveConcentration = (double *) malloc(N * sizeof(double));
	double *belowConcentration = (double *) malloc(N * sizeof(double));
	int *aboveBacilli = (int *) malloc(N * sizeof(int));
	int *belowBacilli = (int *) malloc(N * sizeof(int));
	while (true) {
		MPI_Recv(&stop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		if (stop == 1) {
			break;
		}
		i++;
		double maxRedError = updateConcentration(rank, size, 0, workload,
			concentration, bacilli,
			aboveConcentration, belowConcentration, aboveBacilli, belowBacilli);
		double maxBlackError = updateConcentration(rank, size, 1, workload,
			concentration, bacilli,
			aboveConcentration, belowConcentration, aboveBacilli, belowBacilli);
		double maxError = max(maxRedError, maxBlackError);
		int rWorkload = workload[rank] - workload[rank - 1];
		double sumCon = sumConcentration(concentration, rWorkload);
		double *probability = (double *) malloc(rWorkload * N * sizeof(double));
		int l, r;
		for (l = 0; l < rWorkload; l++) {
			for (m = 0; m < N; m++) {
				probability[N * l + m] = isCandidate() * pow(concentration[N * l + m], eta) / sumCon;
				bacilli[N * l + m] = (rand() % ran) / (ran - 1) < probability[N * l + m]
				? 1 : 0;
			}
		}
		MPI_Isend(concentration, N * rWorkload, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &reportTaskRequest);
		MPI_Isend(bacilli, N * rWorkload, MPI_INT, 0, i, MPI_COMM_WORLD, &reportTaskRequest);
		MPI_Isend(&maxError, 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &reportTaskRequest);
	}
	free(aboveConcentration);
	free(belowConcentration);
	free(aboveBacilli);
	free(belowBacilli);
}

void main(int argc, char **argv) {

	double *concentration;
	int *bacilli;

	MPI_Init(&argc, &argv);

	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		concentration = (double *) malloc (N * N * sizeof(double));
		bacilli = (int *) malloc (N * N * sizeof(int));
		initialize(rank, size, concentration, bacilli);
		int *workload = devideTask(rank, size, concentration, bacilli);
	} else {
		int *workload;
		double *concentration;
		int *bacilli;
	}
	int iter;
	for (iter = 0; iter < numIteration; it++) {
		if (ranl == 0) {
			collectResult(rank, size, concentration, bacilli, workload);
		} else {
			recvTask(rank, size, concentration, bacilli, &workload);
			solveTask(rank, size, concentration, bacilli, workload);	
		}
	}

	MPI_Finalize();

	return;
}