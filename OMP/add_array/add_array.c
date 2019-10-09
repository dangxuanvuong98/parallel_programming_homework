#include <omp.h>
#include <stdio.h>
#include <malloc.h>

void main() {
	int n, *A, *B, *C;
	int i;
	int nthreads;
	int id;
	int ss;

	scanf("%d", &n);
	A = (int *) malloc(n * sizeof(int));
	B = (int *) malloc(n * sizeof(int));
	C = (int *) malloc(n * sizeof(int));

	for (i = 0; i < n; i++) {
		scanf("%d", A + i);
	}

	for (i = 0; i < n; i++) {
		scanf("%d", B + i);
	}

	#pragma omp parallel
		nthreads = omp_get_num_threads();

	ss = n / nthreads;
	printf("nthreads %d\n", nthreads);
	printf("segment size %d\n", ss);

	#pragma omp parallel private(i, id)
	{
		id = omp_get_thread_num();
		for (i = ss*id; i < ss*(id+1); i++) {
			C[i] = A[i] + B[i];
			printf("Thread %d, C[%d] = %d\n", id, i, C[i]);
		}
	}

	for (i = 0; i < n; i++) {
		printf("%d ", C[i]);
	}
}
