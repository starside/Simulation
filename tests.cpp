/*
 * tests.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */

#include <iostream>
#include <math.h>
#include <random>
#include <stdlib.h>
#include "tests.h"
#include "branchedChain.h"

void testCenterOfMass(){
	int numPoints = 14;
	int x = 2;
	branchedChain mol1(numPoints);
	for(int k = 0; k < numPoints; k++) {
		mol1.monomers[k].r.vector[0] = pow(x, (double) k);
		mol1.monomers[k].r.vector[1] = 0;
		mol1.monomers[k].r.vector[2] = 0;
	}
	jVector rcm;
	mol1.findCm(&rcm);
	for(int i = 0; i < 3; i++) {
		std::cout << rcm.vector[i] << std::endl;
	}
	double numerator = 1.0 - pow(x,numPoints);
	double denom = 1.0 - (double)x;
	std::cout << numerator/denom/(double)numPoints << std::endl; //exact via geometrix sum
}

void testRng(){
	stateGen gen,gen2; gen.seed(10);
	std::uniform_real_distribution<double> d1(0.0,1.0);
	for(int i = 0; i < 10; i++){
		std::cout << d1(gen) << std::endl;
	}
	gen2.seed(10);
	gen2.discard(5);
	std::cout << "##############\n";
	for(int i = 0; i < 10; i++){
		std::cout << d1(gen2) << std::endl;
	}
}

void testStructureFactor() {
	stateGen gen;
	branchedChain mol1(CHAINLENGTH);
	mol1.createDendrimerR(CHAINLENGTH,3,5,&gen);
	mol1.sigma = MOL_SIGMA;
	mol1.epsilon = MOL_EPSILON;
	mol1.spacing = MOL_SPACING;
	//mol1.runMC(1000,&gen);

	double rmax = 10.0;
	double qmax = 1.0/0.1; //close range interaction cutoff
	double pointsPerUnit = 50.0;
	size_t sfbins = floor(qmax*pointsPerUnit);
	std::cerr << "Bins is " << sfbins << std::endl;
	double *res = new double[sfbins];

	int numf = 1000;
	for(int i = 0; i < numf; i++) {
		mol1.runMC(500,&gen);
		mol1.findStructureFactorRadial(res,sfbins,qmax);
	}
	for(int i = 0; i < sfbins; i++){
		//q = i*qWidth + qw2; //Use middle of bin
		std::cout << i << ", " << res[i]/(double)numf << std::endl;
	}
	delete[] res;
	exit(0);

}

void testSaveState() {
	int chainlen = 100; int ms = 666;
	int batch = 666;
	dataFP fp;
	stateGen gen; gen.seed(ms);
	branchedChain mol1(chainlen), mol2(chainlen);
	mol1.setLinearPositions(2.0);
	fp = fopen("foo.dat","wb");
	mol1.writeHeader(fp,&gen);
	mol1.writeParameterLine(fp);
	mol1.writeDataLine(fp,&gen);
	fclose(fp);
	fp = fopen("foo.dat","rb");
	mol2.readHeader(fp,&gen);
	mol2.readLine(fp,&gen);
	mol2.readLine(fp,&gen);
	fclose(fp);
	for(int i = 0; i < mol1.numMonomers; i++){
		jVector res;
		vectorSubP(&res, &mol1.monomers[i].r, &mol2.monomers[i].r);
		for(int k = 0; k < 3; k++) {
			std::cout << res.vector[k];
		}
	}
	std::cout << std::endl;
	//Test save and load state
	gen.seed(ms);
	branchedChain mol3(chainlen);
	mol3.createLinear(chainlen);
	mol3.setLinearPositions(1.0);
	mol3.runMC(1000,&gen);
	std::cout << "Rg is " << mol3.findRg() <<  " Batch numer is " << batch << std::endl;
	mol3.saveState("test.dat", &gen, &batch);
	mol3.runMC(1000,&gen);
	std::cout << "Rg is " << mol3.findRg() << " Batch numer is " << batch << std::endl;
	std::cout << "Loading old state " << std::endl;
	mol3.loadState("test.dat",&gen, &batch);
	std::cout << "Rg is " << mol3.findRg() << " Batch numer is " << batch << std::endl;
	mol3.runMC(1000,&gen);
	std::cout << "Rg is " << mol3.findRg() << " Batch numer is " << batch << std::endl;

}

#ifdef USE_GUI

vtkSmartPointer<vtkPolyData> pointsOnSphere(const int nLongitude, const int nAzimuth) {
	const double nL = (double) nLongitude;
	const double nA = (double) nAzimuth;
	//divide up angles based on number of points specified
	const double deltaL = M_PI/(nL+1.0); //the 1.0 prevents plotting at poles.  Imagin nL=1
	const double deltaA = 2.0*M_PI/nA;
	//Find number of data points
	const int numPoints = nLongitude*nAzimuth;
	//Dynamically allocate array
	vtkIdType *pid = new vtkIdType[numPoints]; 	//point ids
	jVector scratch = {{0,0,0}};			//a scratch pad

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	//calculate vectors and fill pid
	for(int l = 0; l < nLongitude; l++){
		for(int a = 0; a < nAzimuth; a++) {
			scratch.vector[0]= sin((l+1.0)*deltaL)*cos(deltaA*a);
			scratch.vector[1]= sin((l+1.0)*deltaL)*sin(deltaA*a);
			scratch.vector[2]= cos((l+1.0)*deltaL);
			pid[l*nAzimuth+a] = points->InsertNextPoint(scratch.vector);
		}
	}
	vertices->InsertNextCell(numPoints,pid);
	vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();
	// Set the points and vertices we created as the geometry and topology of the polydata
  	point->SetPoints(points);
  	point->SetVerts(vertices);
	delete[] pid;
	return point;
}

vtkSmartPointer<vtkPolyData> showMonomers(const branchedChain *mol) {
	const int numPoints = mol->numMonomers;
	vtkIdType *pid = new vtkIdType[numPoints]; 	//point ids
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	//calculate vectors and fill pid
	for(int l = 0; l < numPoints; l++){
		pid[l] = points->InsertNextPoint(mol->monomers[l].r.vector);
	}
	vertices->InsertNextCell(numPoints,pid);
	vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();
	// Set the points and vertices we created as the geometry and topology of the polydata
	point->SetPoints(points);
	point->SetVerts(vertices);
	//point->SetLines(vertices);
	delete[] pid;
	return point;
}

vtkSmartPointer<vtkPolyData> randomPointsOnSphere(const int numPoints) {
	//set up random numbers
	std::default_random_engine generator;
  	std::uniform_real_distribution<double> azDistribution(0,M_PI*2.0);
	std::uniform_real_distribution<double> thetaDistribution(-1.0,1.0);
	//Dynamically allocate array
	vtkIdType *pid = new vtkIdType[numPoints]; 	//point ids
	jVector scratch = {{0,0,0}};			//a scratch pad
	const jVector fv = {{0.0,1.0,0.0}};

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	//calculate vectors and fill pid
	for(int l = 0; l < numPoints; l++){
		double r1 = azDistribution(generator); 		//the azimuthal rotations
		double r3 = azDistribution(generator);
		double r2 = acos(thetaDistribution(generator)); //the "theta" angle

		//scratch.vector[0]= sin(r2)*cos(r1);
		//scratch.vector[1]= sin(r2)*sin(r1);
		//scratch.vector[2]= cos(r2);
		jMatrix randomRot = eulerAngleMatrix(r1,r2,r3); //find the matrix
		const jMatrix * jmp = &randomRot;
		const jVector *fvp = &fv;
		transformVectorP(&scratch, jmp, fvp );	//rotate some vector
		pid[l] = points->InsertNextPoint(scratch.vector);
	}
	vertices->InsertNextCell(numPoints,pid);
	vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();
	// Set the points and vertices we created as the geometry and topology of the polydata
  	point->SetPoints(points);
  	point->SetVerts(vertices);
	delete[] pid;
	return point;
}
#endif