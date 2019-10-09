#include <stdio.h>
#include <malloc.h>

void main() {
  int n, *A, *B, *C;
  int i;

  scanf("%d", &n);
  A = (int *) malloc(n * sizeof(int));
  B = (int *) malloc(n * sizeof(int));
  C = (int *) malloc(n * sizeof(int));

  for (i = 0; i < n; i++) {
    scanf("%d", A + i);
  }

  for (i = 0; i < n; i++) {
    scanf("%d", B + i);
    C[i] = A[i] + B[i];
  }

  for (i = 0; i < n; i++) {
    printf("%d ", C[i]);
  }
}
