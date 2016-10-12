/*
 * data.h
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */

#ifndef DATA_H_
#define DATA_H_

extern "C" {
#include <cblas.h>
}
#include <math.h>
#include <iostream>
#include "definitions.h"
#include "dtypes.h"

#ifndef M_PI
#define M_PI 3.1415
#endif

//higher level function
void basicHistogram(double value, unsigned int *bins, unsigned int nbins, double leftEdge, double rightEdge);
void basicHistogramLabels(double *labels, unsigned int *bins, unsigned int nbins, double leftEdge, double rightEdge);


//Function Declarations
void printjVector(const jVector *vec);
void printjMatrix(const jMatrix *mat);

inline void vectorAddP(jVector *res, const jVector *a, const jVector *b);
inline void vectorSubP(jVector *res, const jVector *a, const jVector *b);
inline double sinc(const double x);

void vectorAddPd(double *res, const double *a, const double *b);
void vectorSubPd(double *res, const double *a, const double *b);
void transformVectorP(jVector* res, const jMatrix *mat, const jVector *vec);
void transformMatrixP(jMatrix *res, const jMatrix *mLeft, const jMatrix *mRight);
jMatrix eulerAngleMatrix(const double r1, const double r2, const double r3);
jMatrix eulerAngleZRotCW(const double theta); //Rotate about z axis
jMatrix eulerAngleXRotCW(const double theta);	//Rotate about x axis
int testSomeRotations(); //some test code to see if matrix and vector transforms work

void vectorSubV(jVector *__restrict__ res, const jVector *__restrict__ a, const jVector *__restrict__ b);

//Add two jVectors
inline void vectorAddP(jVector *res, const jVector *a, const jVector *b) {
	for(int i = 0; i < 3; i++) {
		res->vector[i] = a->vector[i] + b->vector[i];
	}
}

inline void vectorAvgP(jVector *res, const jVector *a, const jVector *b);
inline void vectorAvgP(jVector *res, const jVector *a, const jVector *b) {
	for(int i = 0; i < 3; i++) {
		res->vector[i] = (a->vector[i] + b->vector[i])/2.0;
	}
}

//dot two jVectors
inline double dotP( const jVector *a, const jVector *b);
inline double dotP( const jVector *a, const jVector *b) {
	double v = 0;
	for(int i = 0; i < 3; i++) {
		v += a->vector[i]*b->vector[i];
	}
	return v;
}

//dot two jVectors
inline void scaleP(jVector *res, const jVector *a, const double s);
inline void scaleP(jVector *res, const jVector *a, const double s) {
	for(int i = 0; i < 3; i++) {
		res->vector[i] = a->vector[i]*s;
	}
}

//Add two jVectors
inline void vectorSubP(jVector *__restrict__ res, const jVector * __restrict__ a, const jVector * __restrict__ b) {
	res->vector[0] = a->vector[0] - b->vector[0];
	res->vector[1] = a->vector[1] - b->vector[1];
	res->vector[2] = a->vector[2] - b->vector[2];
}


inline double diffSquared(const jVector *a, const jVector *b);

inline double diffSquared(const jVector *a, const jVector *b) {
	jVector res = {{0,0,0}};
	vectorSubV(&res,a,b);
	return res.vector[0]*res.vector[0] + res.vector[1]*res.vector[1] + res.vector[2]*res.vector[2];
}


//This required 16 bye aligned memory!
 inline void vectorSubV(jVector *__restrict__ res, const jVector *__restrict__ a, const jVector *__restrict__ b) {
	for(int i = 0; i < 3; i++) {
		res->vector[i] = a->vector[i] - b->vector[i];
	}
}

inline double sinc(const double x) {
	if(x <= 0.01) {  //this choice is good to one part in 10^11
		return 1 - x*x/6.0;
	} else {
		return sin(x)/x;
	}
}

#endif /* DATA_H_ */
