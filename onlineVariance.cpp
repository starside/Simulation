#include "onlineVariance.h"

/*
//Headers for testing
#include <random>
#include <iostream>
*/

OnlineVar::OnlineVar() {
	n = 0;
	mean = 0.0;
	m2 = 0.0;
}

void OnlineVar::addValue(double x) {
	double delta = x - mean;
	n += 1;
	mean += delta/(double)n;
	m2 += delta*(x-mean);
}

double OnlineVar::Variance() {
	if(n < 2) {
		return strtod("NAN",NULL);
	} else {
		return m2/(double)(n-1);
	}
}

void OnlineVar::setState(unsigned long int n, double mean, double m2){
	this->n = n;
	this->mean = mean;
	this->m2 = m2;
}

OnlineVar *OLVCombine(OnlineVar *a, OnlineVar*b) {
	OnlineVar * res = new OnlineVar; //allocate new object
	double delta = b->mean - a->mean;
	res->mean = ((double)a->n*a->mean + (double)b->n*b->mean)/(double)(a->n+b->n);
	res->m2 = a->m2 + b->m2 + delta*delta*a->n*b->n/(double)(a->n + b->n);
	res->n = a->n + b->n;
	return res;
}

/*Simple test of the code
int main() {
	std::mt19937_64 gen;
	std::normal_distribution<double> distribution (0.0,1.0);

	OnlineVar l1;
	OnlineVar l2;
	OnlineVar lT;

	for(int i =0; i < 50000; i++) {
		double nv1 = distribution(gen);
		double nv2 = distribution(gen);
		double nv3 = distribution(gen);
		l1.addValue(nv1);
		l2.addValue(nv2);l2.addValue(nv3);
		lT.addValue(nv1); lT.addValue(nv2); lT.addValue(nv3);
	}
	std::cout << "l1 has Variance " << l1.Variance() <<" l2 has variance " << l2.Variance() << "  lT has variance " << lT.Variance() << std::endl;
	std::cout << "l1 has mean " << l1.mean <<" l2 has mean " << l2.mean << "  lT has mean " << lT.mean << std::endl;
	std::cout << "l1 has m2 " << l1.m2 <<" l2 has m2 " << l2.m2 << "  lT has m2 " << lT.m2 << std::endl;

	OnlineVar *olc = OLVCombine(&l1,&l2);
	std::cout << "Combined l1 and l2 has variance " << olc->Variance() << std::endl;
	std::cout << "Combined l1 and l2 has mean " << olc->mean << std::endl;
	std::cout << "Combined l1 and l2 has m2 " << olc->m2 << std::endl;
	return 0;
}
*/