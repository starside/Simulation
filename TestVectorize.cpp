#include <iostream>
#include <cstdint>
#define MAXBRANCHES 1000000000
typedef struct {double matrix[9];} jMatrix;
typedef struct {double vector[4];} jVector;
typedef struct {jVector r; int branch[MAXBRANCHES]; int zeroBranch; } jMonomer;

 void vectorSubV(jVector *__restrict__ res, const jVector *__restrict__ a, const jVector *__restrict__ b);


 void vectorSubV(jVector *__restrict__ res, const jVector *__restrict__ a, const jVector *__restrict__ b) {
	double *x = (double *)__builtin_assume_aligned(a->vector, 16);
	double *y = (double *)__builtin_assume_aligned(b->vector, 16);
	double *z = (double *)__builtin_assume_aligned(res->vector, 16);
	for(int i = 0; i < 3; i++) {
		z[i] = x[i] - y[i];
		z[i] = z[i]*z[i];
	}
}

int main() {
	unsigned int a,b;
	jVector *mj = new jVector[3];
	std::cout << (uint64_t)mj % 16 << std::endl;
	a = clock();
	for(uint64_t i = 0; i < MAXBRANCHES; i++){
		vectorSubV(&mj[0], &mj[1], &mj[2]);
	}
	b = clock();
	std::cout << b-a << std::endl;
	return 0;
}