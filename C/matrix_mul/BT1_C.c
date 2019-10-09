#include <stdio.h>
#include <malloc.h>

void main() {
	int m, n, p;
	int **A, **B, **C;
	int i, j, k;

	scanf("%d%d%d", &m, &n, &p);

	A = (int **) malloc(m * sizeof(int *));
	B = (int **) malloc(n * sizeof(int *));
	C = (int **) malloc(m * sizeof(int *));

	for (i = 0; i < m; i++) {
		A[i] = (int *) malloc(n * sizeof(int));
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

	for (i = 0; i < m; i++) {
		C[i] = (int *) malloc(p * sizeof(int));
		for (j = 0; j < p; j++) {
			C[i][j] = 0;
			for (k = 0; k < n; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
			printf("%d ", C[i][j]);
		}
		printf("\n");
	}
}