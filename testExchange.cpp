#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <Eigen/Geometry>
#include <cmath>

double circleArcSin(double sine, double cosine);
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
Eigen::Vector3d monomerSub(Eigen::Matrix3Xd monomers, const int a, const int b);
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

int calculatePerm(Eigen::Matrix3Xd in, int frame[3]);
int calculatePerm(Eigen::Matrix3Xd in, int frame[3]){
	//int frame[3] = {1,2,3};
	int basei = 0; // base monomer
	int basej = frame[0]; // reference monomer
	int numMonomers = 3;
	// A list of pairs that contains a monomer and the angle on the unit circle
	std::vector<std::pair<int, double>> perm;
	//Calculate cross products
	Eigen::Vector3d cp, vec1,vec2;
	for(int i = 0; i < numMonomers; i++) {
			cp = monomerSub(in, basei, basej).cross(monomerSub(in,basei,frame[i]));
			vec1 =  monomerSub(in, basei, basej);
			vec2 = monomerSub(in,basei,frame[i]);
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
	return perm_number;
}

int main() {
	int numMonomers = 4;
	Eigen::Matrix3Xd in(3, numMonomers);
	for(int row = 0; row < 3; row++){
		for(int col=0; col < numMonomers; col++){
			in(row,col) = 0;
		}
	}
	double thetas[3] = {0, 0.5, 1.2};
	for(int m = 0; m < 3; m++){
		in(0, 1 + m) = std::cos(thetas[m]+M_PI/2); // x
		in(1, 1 + m) = std::sin(thetas[m]+M_PI/2); // y
	}

	int frame[3] = {1,2,3};
	int perm_number = calculatePerm(in, frame);
	std::cout << "Permuation is " << perm_number << std::endl;


	return 0;
}
