#include <omp.h>
#include <stdio.h>

void main() {
	omp_set_num_threads(10); //Lenh nay dung truoc #pragma
	#pragma omp parallel
		printf("Hello from thread %d, nthreads %d\n",
			omp_get_thread_num(),
			omp_get_num_threads());
}