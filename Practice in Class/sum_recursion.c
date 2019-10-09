#include <omp.h>
#include <stdio.h>
#include <malloc.h>

int cal_sum(int *a, int l, int r) {
	if (l == r) {
		printf("%d %d %d\n", l, r, a[l]);
		return a[l];
	}
	int m = (l+r) / 2;
	int s, s0, s1;
	
	#pragma omp parallel
	#pragma omp single nowait
	{
		#pragma task shared(s0)
		s0 = cal_sum(a, l, m);
		#pragma task shared(s1)
		s1 = cal_sum(a, m+1, r);
		#pragma task wait
		s = s0 + s1;
	}

	printf("%d %d %d\n", l, r, s);
	return s;
}

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