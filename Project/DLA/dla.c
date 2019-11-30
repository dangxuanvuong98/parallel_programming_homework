#include "mpi.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>

#define TAG_ABOVE_BOUNDARY 1
#define TAG_BELOW_BOUNDARY 2
#define TAG_STOP 3
#define TAG_CONCENTRATION 4
#define TAG_BACILLUS 5
#define TAG_WORKLOAD 6
#define TAG_ERROR 7
#define TAG_SUM_CONCETRATION 8

const int numIteration = 1000;
const double epsilon = 0.001;
const double omega = 1.5;
const double eta = 2.0;
const int N = 40;
const int[] dl = {1, 0, -1, 0};
const int[] dm = {0, 1, 0, -1};

double max(double a, double b) {
	return a > b ? a : b;
}

double _SOR(double *concentration, int l, int m) {
	return (omega / 4) * (
		concentration[N * (l - 1) + m] +
		concentration[N * (l + 1) + m] +
		concentration[N * l + ((m + 1 + N) % N)] +
		concentration[N * l + ((m - 1 + N) % N)]) +
		(1 - omega) * concentration[N * l + m];
}

void log(int iter, double *concentration, int *bacillus) {
	char *logFile = (char *) malloc (1024 * sizeof(char));
	sprintf("log_N=%d_eta=%llf_iter=%d", N, eta, iter);
	FILE *f = fopen(logFile, "w");
	free(logFile);
	int l, m;
	for (l = 0; l < N; l++) {
		for (m = 0; m < N; m++) {
			fprintf(f, "%llf ", concentration[N * l + m]);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	for (l = 0; l < N; l++) {
		for (m = 0; m < N; m++) {
			fprintf(f, "%d ", bacillus[N * l + m]);
		}
		fprintf(f, "\n");
	}
}


void boundariesConcentrationExchange(int rank, int size,
	double *concentration, int *bacillus, int *workload) {
	
	MPI_Status status;

	if (rank > 1) {
		MPI_Send(concentration + N, N, MPI_DOUBLE,
			rank - 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD);
	}
	if (rank < size - 1) {
		MPI_Send(concentration + N * workload, N, MPI_DOUBLE,
			rank + 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD);
	}
	if (rank > 1) {
		MPI_Recv(concentration + N * (workload + 1), N, MPI_DOUBLE,
			rank - 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD, &status);
	}
	if (rank < size - 1) {
		MPI_Recv(concentration, N, MPI_DOUBLE,
			rank + 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD, &status);
	}
}

void boundariesBacillusExchange(int rank, int size,
	double *concentration, int *bacillus, int *workload) {
	if (rank > 1) {
		MPI_Send(bacillus + N, N, MPI_INT,
			rank - 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD);
	}
	if (rank < size - 1) {
		MPI_Send(bacillus + N * workload, N, MPI_INT,
			rank + 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD);
	}
	if (rank > 1) {
		MPI_Recv(bacillus + N * (workload + 1), N, MPI_INT,
			rank - 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD, &status);
	}
	if (rank < size - 1) {
		MPI_Recv(bacillus, N, MPI_INT,
			rank + 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD, &status);
	}
}

double updateConcentration(int rank, int size,
	double *concentration, int *bacillus, int *workload, int rb) {

	double maxError = 0;
	boundariesConcentrationExchange(rank, size, concentration, bacillus, workload);

	int l, m;
	int rWorkload = workload[rank] - workload[rank - 1];
	for (l = 1; l <= rWorkload; l++) {
		for (m = 0; m < N; m++) {
			int tmp = (l - 1) + r + N * workload[rank - 1];
			if (tmp % 2 == rb) {
				double oldConcentration = concentration[N * l + m];
				concentration[N * l + m] = _SOR(concentration, l, m);
				if (bacillus[N * l + m]) {
					concentration[N * l + m] = 0;
				}
				if (rank == 1 && l == 1) {
					concentration[N * l + m] = 0;
				}
				if (rank == size - 1 && l == rWorkload) {
					concentration[N * l + m] = 1;
				}
				maxError = max(maxError, fabs(concentration[N * l + m] - oldConcentration));
			}
		}
	}
	return maxError;
}

void initialize(int rank, int size, double **concentration, int **bacillus, int **workload) {

	*concentration = (double *) malloc (N * N * sizeof(double));
	*bacillus = (int *) malloc (N * N * sizeof(int));
	*workload = (int *) malloc (size * sizeof(int));

	int l, m;
	for (m = 0; m < N; m++) {
		for (l = 0; l < N - 1; l++) {
			(*concentration)[N * l + m] = 0;
			(*bacillus)[N * l + m] = 0;
		}
		(*concentration)[N * (N - 1) + m] = 1;
		(*bacillus)[N * (N - 1) + m] = 0;
	}
	(*bacillus)[N / 2] = 1;

	int avg_workload = N / (size - 1);
	int over_workload = N % (size - 1);

	(*workload)[0] = 0;
	for (r = 1; r < size; r++) {
		(*workload)[r] = (*workload)[r - 1] + avg_workload;
		if (r <= over_workload) {
			(*workload)[r] += 1;
		}
	}
}

void divideTask(int rank, int size, double *concentration, int *bacillus, int *workload) {
	int stop = 0;
	int r;
	for (r = 1; r < size; r++) {
		int rWorkload = workload[r] - workload[r - 1];
		MPI_Send(workload, size, MPI_INT,
			r, TAG_WORKLOAD, MPI_COMM_WORLD);
		MPI_Send(concentration + N * workload[r - 1], N * rWorkload, MPI_DOUBLE,
			r, TAG_CONCENTRATION, MPI_COMM_WORLD);
		MPI_Send(bacillus + N * workload[r - 1], N * rWorkload, MPI_INT,
			r, TAG_BACILLUS, MPI_COMM_WORLD);
		MPI_Send(&stop, 1, MPI_INT,
			r, TAG_STOP, MPI_COMM_WORLD);
	}
}

void recvTask(int rank, int size, double **concentration, int **bacillus, int **workload) {

	MPI_Status status;

	*workload = (int *) malloc (size * sizeof(int));
	MPI_Recv(*workload, size, MPI_INT,
		0, TAG_WORKLOAD, MPI_COMM_WORLD, &status);
	int rWorkload = (*workload)[rank] - (*workload)[rank - 1];

	*concentration = (double *) malloc (N * (rWorkload + 2) * sizeof(double));
	MPI_Recv(*concentration + N, N * rWorkload, MPI_DOUBLE,
		0, TAG_CONCENTRATION, MPI_COMM_WORLD, &status);

	*bacillus = (int *) malloc (N * (rWorkload + 2) * sizeof(int));
	MPI_Recv(*bacillus + N, N * rWorkload, MPI_INT,
		0, TAG_BACILLUS, MPI_COMM_WORLD, &status);

	if (rank == 1) {
		int m;
		for (m = 0; m < N; m++) {
			(*bacillus)[m] = 0;
		}
	} else if (rank == size - 1) {
		int m;
		for (m = 0; m < N; m++) {
			(*bacillus)[N * (workload + 1) + m] = 0;
		}
	}
}

void collectConcentration(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	double error, maxError;
	int stop = 0;
	int i = 0;
	MPI_Status status;

	int r;
	while (true) {
		maxError = 0;
		i++;
		for (r = 1; r < size; r++) {
			int rWorkload = workload[r] - workload[r - 1];
			MPI_Recv(concentration, N * rWorkload, MPI_DOUBLE,
				r, TAG_CONCENTRATION, MPI_COMM_WORLD, &status);
			MPI_Recv(&error, 1, MPI_DOUBLE,
				r, TAG_ERROR, MPI_COMM_WORLD, &status);
			maxError = max(maxError, error);
		}
		if (maxError < epsilon) {
			stop = 1;
		}
		for (r = 1; r < size; r++) {
			MPI_Send(&stop, 1, MPI_DOUBLE,
				r, TAG_STOP, MPI_COMM_WORLD);
		}
	}
}

void solveDiffusionEquation(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	int stop;
	int i = 0;
	MPI_Status status;

	while (true) {
		MPI_Recv(&stop, 1, MPI_INT,
			0, TAG_STOP, MPI_COMM_WORLD, &status);
		if (stop == 1) {
			break;
		}
		i++;

		double maxRedError = updateConcentration(rank, size,
			concentration, bacillus, workload, 0);
		double maxBlackError = updateConcentration(rank, size,
			concentration, bacillus, workload, 1);

		double maxError = max(maxRedError, maxBlackError);

		MPI_Send(concentration, N * rWorkload, MPI_DOUBLE,
			0, TAG_CONCENTRATION, MPI_COMM_WORLD);
		MPI_Send(&maxError, 1, MPI_DOUBLE,
			0, TAG_ERROR, MPI_COMM_WORLD);
	}
}

void calculateGlobalCandidateConcentration(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	double localCandidateConcentration;
	double globalCandidateConcentration = 0;

	int r;
	for (r = 1; r < size; r++) {
		MPI_Recv(&localCandidateConcentration, 1, MPI_DOUBLE,
			r, TAG_SUM_CONCETRATION, MPI_COMM_WORLD, &status);
		globalCandidateConcentration += localCandidateConcentration;
	}
	
	for (r = 1; r < size; r++) {
		MPI_Send(&globalCandidateConcentration, 1, MPI_DOUBLE,
			r, TAG_SUM_CONCETRATION, MPI_COMM_WORLD);
	}
}

void randomGrow(int rank, int size,
	double *concentration, double *bacillus, int *workload) {

	boundariesBacillusExchange(rank, size,
		concentration, bacillus, workload);

	double localCandidateConcentration = 0;
	double globalCandidateConcentration;

	int l, m, i, tmp;
	for (l = 1; l <= workload; l++) {
		for (m = 0; m < N; m++) {
			tmp = 0;
			for (i = 0; i < 4; i++) {
				tmp = max(tmp, bacillus[N * (l + dl[i]) + (m + dm[i])]);
			}
			if (bacillus[N * l + m] == 1) {
				tmp = 0;
			}
			localCandidateConcentration += tmp * concentration[N * l + m];
			bacillus[N * l + m] = -tmp;
		}
	}

	MPI_Send(&localCandidateConcentration, 1, MPI_DOUBLE,
		0, TAG_SUM_CONCETRATION, MPI_COMM_WORLD);
	MPI_Recv(&globalCandidateConcentration, 1, MPI_DOUBLE,
		0, TAG_SUM_CONCETRATION, MPI_COMM_WORLD);

	double growProbability;
	double grow;
	for (l = 1; l <= workload; l++) {
		for (m = 0; m < N; m++) {
			tmp = -bacillus[N * l + m];
			bacillus[N * l + m] = 0;

			growProbability = tmp * concentration[N * l + m] / globalCandidateConcentration;
			grow = 1.0 * (rand() % (N * N)) / (N * N);
			if (grow < growProbability) {
				bacillus[N * l + m] = 1;
				concentration[N * l + m] = 0;
			}
		}
	}
}

void bacillusGrow(int rank, int size, double *concentration, int *bacillus, int *workload) {

	int rWorkload = workload[rank] - workload[rank - 1];

	if (rank == 0) {
		calculateGlobalCandidateConcentration(rank, size, concentration, bacillus, workload);
	} else {
		randomGrow(rank, size, concentration, bacillus, workload);
	}
}

void diffusionSteady(int rank, int size, double *concentration, int *bacillus,
	int *workload, int iter) {

	if (rank == 0) {
		collectConcentration(rank, size, concentration, bacillus, workload);
	} else {
		solveDiffusionEquation(rank, size, concentration, bacillus, workload);
	}

	if (rank == 0) {
		log(iter, concentration, bacillus);
	}
}

void main(int argc, char **argv) {

	MPI_Init(&argc, &argv);

	srand(time(0));

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double *concentration;
	int *bacillus;
	int *workload;

	if (rank == 0) {
		initialize(rank, size, &concentration, &bacillus, &workload);
		divideTask(rank, size, concentration, bacillus, workload);
	} else {
		recvTask(rank, size, &concentration, &bacillus, &workload);
	}
	diffusionSteady(rank, size, concentration, bacillus, workload, 0);

	int iter;
	for (iter = 1; iter <= numIteration; iter++) {
		bacillusGrow(rank, size, concentration, bacillus, workload);
		diffusitonSteady(rank, size, concentration, bacillus, workload, iter);
	}

	free(concentration);
	free(bacillus);
	free(workload);

	MPI_Finalize();

	return;
}