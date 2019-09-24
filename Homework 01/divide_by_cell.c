#include <omp.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>

void main() {
	int m, n, p;
	int **A, **B, **C;
	int i, j, k;
	int *bn, *bp, *b;
	int ss;
	int nthreads;
	int id;

	// Read input and alloc
	scanf("%d%d%d", &m, &n, &p);
	A = (int **) malloc(m * sizeof(int *));
	B = (int **) malloc(n * sizeof(int *));
	C = (int **) malloc(m * sizeof(int *));
	for (i = 0; i < m; i++) {
		A[i] = (int *) malloc(n * sizeof(int));
		C[i] = (int *) malloc(k * sizeof(int));
		for (j = 0; j < n; j++) {
			scanf("%d", &A[i][j]);
		}
	}
	for (i = 0; i < n; i++) {
		B[i] = (int *) malloc(p * sizeof(int));
		for (j = 0; j < p; j++) {
			scanf("%d", &B[i][j]);
		}
	}

	// Divide
	#pragma omp parallel
		nthreads = omp_get_num_threads();
	b = (int*) malloc((nthreads+1) * sizeof(int));
	bn = (int*) malloc((nthreads+1) * sizeof(int));
	bp = (int*) malloc((nthreads+1) * sizeof(int));
	ss = sqrt(nthreads);
	if (ss * ss < nthreads) {
		ss++;
	}
	for (i = 0; i < ss; i++) {
		
	}

	// Conquer
	#pragma omp parallel private(i, j, k, id)
	{
		id = omp_get_thread_num();
		for (i = 0; i < m; i++) {
			for (j = b[id]; j < b[id+1]; j++) {
				C[i][j] = 0;
				for (k = 0; k < n; k++) {
					C[i][j] += A[i][k] * B[k][j];
				}
				printf("Thread %d calculate element (%d, %d)\n", id, i, j);
			}
		}
	}

	for (i = 0; i < m; i++) {
		for (j = 0; j < p; j++) {
			printf("%d ", C[i][j]);
			}
		printf("\n");
	}
}