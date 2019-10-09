#include <omp.h>
#include <stdio.h>

void main() {
	omp_set_num_threads(100); //Lenh nay dung truoc #pragma
	#pragma omp parallel
	{
		int id, x;
		id = omp_get_thread_num();
		x = 10*id;
		printf("\n");
		printf("Hello from thread %d,x = %d", id, x);
		printf("\n");
	}
}