#ifndef ONLINE_VAR_H
#define ONLINE_VAR_H

#include <math.h>
#include <stdlib.h>

class OnlineVar {
	public:
		unsigned long int n;
		double mean;
		double m2;
		OnlineVar();
		void addValue(double x);
		double Variance();
		void setState(unsigned long int n, double mean, double m2);
};

OnlineVar *OLVCombine(OnlineVar *a, OnlineVar*b); //combines two Objects, and returns a third
#endif