#include <omp.h>
#include <stdio.h>
#include <malloc.h>

void main() {
	int n, *A;
	int i;

	scanf("%d", &n);
	A = (int *) malloc(n * sizeof(int));

	for (i = 0; i < n; i++) {
		scanf("%d", A + i);
	}

	int sum = cal_sum(A, 0, n-1);

	printf("%d\n", sum);
}