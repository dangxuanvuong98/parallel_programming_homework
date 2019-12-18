#include <stdio.h>

__global__ void mykernel(void) {
	printf("Hello World!\n");
}

int main(void) {
	mykernel<<<1, 1>>>();
	printf("Hello World!\n");
	return 0;
}