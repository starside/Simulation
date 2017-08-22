/*
 * data.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */
#include "data.h"
#include <cstring>
#include <math.h>


//Function Definitions
//calculates a very simple histogram.  Not efficient at all
void basicHistogram(double value, unsigned int *bins, unsigned int nbins, double leftEdge, double rightEdge){
	double nv = (value - leftEdge)/(rightEdge - leftEdge);
	if(nv > 1.0 or nv < 0.0) {
		//std::cerr << "Hitogram value out of bounds" << std::endl;
		//exit(0);
		return;
	}
	double dn = 1.0/(double) nbins;
	unsigned int i = floor(nv/dn);
	//add value to histogram
	bins[i] += 1;
}

void basicHistogramLabels(double *labels, unsigned int *bins, unsigned int nbins, double leftEdge, double rightEdge){
	double dn = (rightEdge - leftEdge)/(double) nbins;
	//calculate labels
	for(int k = 0; k < nbins; k++) {
		labels[k] = dn/2.0 + (double)k*dn + leftEdge;
		bins[k] = 0;
	}
}



void printjVector(const jVector *vec){
	for(int i=0; i < 3; i++) {
		std::cout << vec->vector[i] << ",";
	}
	std::cout << "\n";
}

//prints a column major order (Fortran Style) matrix
void printjMatrix(const jMatrix *mat){
	int index;
	for(int row = 0; row < 3; row++){
		for(int col = 0; col < 3; col++){
			index = col*3 + row;
			std::cout << mat->matrix[index] << ",";
		}
		std::cout << "\n";
	}
}

void vectorAddPd(double *res, const double *a, const double *b, const double bsign) {
	res[0] = a[0] + bsign*b[0];
	res[1] = a[1] + bsign*b[1];
	res[2] = a[2] + bsign*b[2];
}

/*The return value is a location in memory.  I plan on rotating many
vectors, so not using the stack seems beneficial*/
void transformVectorP(jVector* res, const jMatrix *mat, const jVector *vec){
	//The LDA paramater is the stride.  In other words, how many
	//items per column.  It is used to map to linear memory
	cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, mat->matrix, 3,
		vec->vector, 1, 0, res->vector, 1);
}

/*The return value is a location in memory.  Perform matrix multiplication*/
void transformMatrixP(jMatrix *res, const jMatrix *mLeft, const jMatrix *mRight) {
	cblas_dgemm(	CblasColMajor, CblasNoTrans,
			CblasNoTrans, 3, 3,
			3, 1.0, mLeft->matrix,
			3, mRight->matrix, 3,
			0, res->matrix, 3);
}

//I am using Jose and Saletan p.528 for this definition
jMatrix eulerAngleZRotCW(const double theta){
	jMatrix result;
	/*
	cos t -sin t 0
	sin t  cos t 0
	0      0     1
	*/
	result.matrix[0] = cos(theta);
	result.matrix[1] = sin(theta);
	result.matrix[2] = 0;
	result.matrix[3] = -sin(theta);
	result.matrix[4] = cos(theta);
	result.matrix[5] = 0;
	result.matrix[6] = 0; result.matrix[7] = 0; result.matrix[8] = 1.0;
	return result;
}

//I am using Jose and Saletan p.528 for this definition
jMatrix eulerAngleXRotCW(const double theta){
	jMatrix result;
	/*
	1      0     0
	0    cos t  -sin t
	0    sin t   cos t
	*/
	result.matrix[0] = 1;result.matrix[1] = 0;result.matrix[2] = 0;
	result.matrix[3] = 0;
	result.matrix[4] = cos(theta);
	result.matrix[5] = sin(theta);
	result.matrix[6] = 0;
	result.matrix[7] = -1.0*sin(theta);
	result.matrix[8] = cos(theta);
	return result;
}

/*This is a body->lab frame transform
Pretty wasteful with the stack, but who cares.  Seems
a little safe than dynamic allocation
x_lab = M x_body.  To go the other direction, swap
psi and phi, then invert all signs.  Needs testing
*/
jMatrix eulerAngleMatrix(const double phi, const double theta, const double psi){
	jMatrix mPhi, mTheta, mPsi,mEuler,mTemp;
	mPhi   = eulerAngleZRotCW(phi);
	mTheta = eulerAngleXRotCW(theta);
	mPsi   = eulerAngleZRotCW(psi);
	transformMatrixP(&mTemp,  &mTheta, &mPsi);
	transformMatrixP(&mEuler, &mPhi,   &mTemp);
	return mEuler;
}


int testSomeRotations(){
	jVector test = { {0,0,1} };
	jVector res = { {6,6,6} };
	jMatrix resM;
	jMatrix rot1 = eulerAngleXRotCW(45*M_PI/180.0);
	jMatrix rot2 = eulerAngleZRotCW(45*M_PI/180.0);
	transformVectorP(&res, &rot1, &test);

	transformMatrixP(&resM,&rot2,&rot1);
	//printjVector(&res);
	std::cout << "----M1----\n";
	printjMatrix(&rot1);
	std::cout << "----M2----\n";
	printjMatrix(&rot2);
	std::cout << "---M1*M2---\n";
	transformMatrixP(&resM,&rot1,&rot2);
	printjMatrix(&resM);
	std::cout << "---M2*M1---\n";
	transformMatrixP(&resM,&rot2,&rot1);
	printjMatrix(&resM);
	return 0;
}

/* The Following block of code is used for finding which permuation the inner
   ring of a tri-functional dendrimer is in. It is NOT general */

double circleArcSin(double sine, double cosine){
	// If[cs >= 0,If[as>0,as, 2\[Pi]+as] ,\[Pi]-as ]
	double as = std::asin(sine);
	if(cosine >= 0) {
		if(as >= 0) {
			return as;
		} else{
			return 2*M_PI + as;
		}
	} else {
		return M_PI - as;
	}
}

//Simple way to subtract 2 monomers
Eigen::Vector3d monomerSub(Eigen::Matrix3Xd monomers, const int a, const int b){
	Eigen::Vector3d res(0,0,0); // 3 column vector
	for(int c = 0; c < 3; c++) {
		res(c) = monomers(c,b) - monomers(c,a);
	}
	return res.normalized();
}

bool sort_function(std::pair<int, double> a, std::pair<int, double> b) {
	return a.second < b.second;
}

int calculatePerm(Eigen::Matrix3Xd in, int frame[]){
	//int frame[3] = {1,2,3};
	int basei = frame[0]; // base monomer
	int basej = frame[1]; // reference monomer
	int numMonomers = 4;
	// A list of pairs that contains a monomer and the angle on the unit circle
	std::vector<std::pair<int, double>> perm;
	//Calculate cross products
	Eigen::Vector3d cp, vec1,vec2;
	for(int i = 1; i < numMonomers; i++) {
			vec1 =  monomerSub(in, basei, basej);
			vec2 = monomerSub(in,basei,frame[i]);
			cp = vec1.cross(vec2);
			double cosa = vec1.dot(vec2);
			double sina = cp(2); //sin of angle
			perm.push_back(std::make_pair(frame[i],circleArcSin(sina, cosa)));
	}

	//Sort the permuation
	std::sort(perm.begin(), perm.end(), sort_function);
	int perm_number = 0;
	//Convert permuation to an integer.
	for(std::vector<std::pair<int,double>>::iterator it = perm.begin(); it != perm.end(); ++it){
		perm_number += it->first;
		perm_number *= 10;
	}
	perm_number /= 10; //divide out extra 10

	if(perm_number == 123){
		return 0;
	} else if (perm_number == 132){
		return 1;
	} else{
		std::cerr << "Got unexpected permuation" << std::endl;
		exit(0);
	}

}
