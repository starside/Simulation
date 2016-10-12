#ifndef DTYPES_H_
#define DTYPES_H_

#include "definitions.h"

typedef struct log_line {
	int N;
	double rog2_sum ; //square radius of gyration sum
	double de_sum; //energy derivative
	double rmax; //maximum extent of density
	double density_sum[DENSITY_BINS];
	double spacing;
	double epsilon;
	double sigma;
	//Data for online variance calculations
	unsigned long int rg2_n; //radius of gyration squared
	double rg2_mean;
	double rg2_m2;
	unsigned long int de_n; //energy derivative
	double de_mean;
	double de_m2;
	int edgeCount;

	unsigned long int ulongcandle;
} LogLine;

//Data Structures
typedef struct {double matrix[9];} jMatrix;
typedef struct {double vector[3]; } jVector; //4th is for padding.  Keep 16 byte alignment
typedef struct {
	jVector r; 
	int branch[MAXBRANCHES]; //the plus 2 is a placeholdher for pointers to phantom monomers 
	int zeroBranch;
	double lastInteraction; 
} jMonomer;
const double zeroVector[3] = {0,0,0};

#endif