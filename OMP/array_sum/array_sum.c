#include <omp.h>
#include <stdio.h>
#include <malloc.h>

void main() {
	int n, *A;
	int nthreads, ss;
	int sum = 0;
	int i;
	scanf("%d", &n);
	A = (int *) malloc(n * sizeof(int));

	for (i = 0; i < n; i++) {
		scanf("%d", A + i);
	}

	#pragma omp parallel
	{
		int i, id;
		nthreads = omp_get_num_threads();
		id = omp_get_thread_num();
		ss = n / nthreads;
		for (i = id*ss; i < (id+1)*ss; i++) {
			sum += A[i];
		}
	}

	printf("%d\n", sum);
}