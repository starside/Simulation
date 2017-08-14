/*
 * branchedChain.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */

#include <stdlib.h>
#include <string.h>
#include <queue>
#include <random>
#include <iomanip>
#include <unistd.h>
#include <sys/types.h>
#include <sstream>
#include <fstream>
#include "branchedChain.h"
#include "definitions.h"

inline double unitPot(const double d2, const double sig);

//This function is too slow to use in the same place as find density!
//Find the structure factor, and store in double array of length rlen
//qmax is the largest q value.
void branchedChain::findStructureFactorRadial(double *result, size_t rlen, double qmax){
	double rij;
	double qWidth = qmax/(double)rlen;
	double qw2 = qWidth/2.0;
	double dnm = 2.0/(double)numMonomers;
	double q;
	for(int i = 0; i < numMonomers; i++) {
		for(int j=i+1; j < numMonomers; j++) {
			rij = sqrt(diffSquared(&monomers[i].r, &monomers[j].r));
			for(int z = 0; z < rlen; z++){ //compute result
				q = z*qWidth + qw2; //Use middle of bin
				result[z] += dnm*sinc(rij*q);
			}
		}
	}
	for(int i = 0; i < rlen; i++){ //add on diagonal terms
		result[i] += 1.0;
	}
}

double branchedChain::findRg() {
	//return false;  //disable SARW
    double rg = 0.0; 
	int nm = this->numMonomers - this->numPhantoms;
	jMonomer *monomers = this->monomers;
	for(int i = 0; i < nm; i++) {
		for(int j=i+1; j < nm; j++){
			double d2 = diffSquared(&monomers[i].r, &monomers[j].r);
			rg += d2;
		}
	}
	return rg/((double)(nm*nm));
}

double branchedChain::findRg(const int dim) {
	//return false;  //disable SARW
    double rg = 0.0; 
    jVector res = {{0,0,0}};
	int nm = this->numMonomers - this->numPhantoms;
	jMonomer *monomers = this->monomers;
	for(int i = 0; i < nm; i++) {
		for(int j=i+1; j < nm; j++){
			vectorSubV(&res,&monomers[i].r, &monomers[j].r);
			double d2 = res.vector[dim]*res.vector[dim];
			rg += d2;
		}
	}
	return rg/((double)(nm*nm));
}

//finds center of mass
void branchedChain::findCm(jVector *res) {
	//return false;  //disable SARW
    jVector rcm = {{0,0,0}}; //zero out the vector
	int nm = this->numMonomers - this->numPhantoms;
	jMonomer *monomers = this->monomers;
	for(int i = 0; i < nm; i++) {
		vectorAddP(&rcm, &rcm,&monomers[i].r);
	}
	double N = (double)nm;
	for(int i = 0; i < 3; i++){ //divide by N
		rcm.vector[i] = rcm.vector[i]/N;
	}
	for(int i = 0; i < 3; i++) { //copy result
		res->vector[i] = rcm.vector[i];
	}
}

//Methods
void branchedChain::printConnections(){
	for(int i = 0; i < numMonomers; i++) {
		std::cout << "Monomer "<< i <<": ";
		for(int b = 0; b < MAXBRANCHES; b++) {
			std::cout << monomers[i].branch[b] << ", ";
		}
		std::cout <<"::"<< monomers[i].zeroBranch << std::endl;
	}
}

void branchedChain::linkTwo(const int src, const int dest) {
	std::cout << "Linking: " << src << " to" << dest << std::endl;
	if(src >= numMonomers || dest >= numMonomers ||
		src < 0 || dest < 0) { std::cerr << "1) Cannot link two monomers "<< src <<","<<dest<<".  You got a bug!\n"; }
	//Find first self-link
	int freeLinkS = -1; int freeLinkD = -1;
	for(int i = 0; i < MAXBRANCHES; i++) {
		if(monomers[src].branch[i] == src) { //linked to self
			freeLinkS = i;
		}
		if(monomers[dest].branch[i] == dest) { //linked to self
			freeLinkD = i;
		}
		if(freeLinkD == -1 && freeLinkS == -1){
			break;
		}

	}
	if(freeLinkS == -1 || freeLinkD == -1) {
		std::cerr << "2) Cannont link two monomers "<< src <<","<<dest<<" FreelinkS and D " << freeLinkS <<"," << freeLinkD << ".  You got a bug!\n";
		exit(0);
	}
	//Link up the monomers
	monomers[src].branch[freeLinkS] = dest; //0 is "forward"
	monomers[dest].branch[freeLinkD] = src; //1 is "backward"
}

void branchedChain::unlinkTwo(const int src, const int dest) {
	if(src >= numMonomers || dest >= numMonomers ||
		src < 0 || dest < 0) { std::cerr << "Cannont unlink two monomers.  You got a bug!\n"; }
	//Find first self-link
	for(int i = 0; i < MAXBRANCHES; i++) {
		if(monomers[src].branch[i] == dest) { //linked to dest
			monomers[src].branch[i] = src;  //unlink
		}
		if(monomers[dest].branch[i] == src) { //linked to self
			monomers[dest].branch[i] = dest;  //unlink
		}
	}
}

void branchedChain::saveTopofile(char *fname) {
	std::ofstream fs;
	fs.open(fname,std::ios::out);
	for(int row = 0; row < numMonomers; row++){
		for(int col = 0; col < numMonomers; col++) {
			int bv = 0;
			for(int b=0; b < MAXBRANCHES; b++){
				if(monomers[row].branch[b] == col && monomers[row].branch[b] != row) {
					bv =1;
				}
			}
			fs << bv;
			if(col < numMonomers - 1) {fs << ",";}
		}
		fs << std::endl;
	}
	fs.close();
}

/*Creates a linear topology.  Does not position monomers*/
void branchedChain::createLinear(int cl) {
	monomers[0].zeroBranch = 0;
	for(int i = 0; i < cl-1; i++) {
		linkTwo(i,i+1);
		monomers[i+1].zeroBranch = i;
	}
	this->buildParentMatrix();
}
//double r[3]
void branchedChain::setLinearPositions(const double spacing) {
	for(int i = 0; i < numMonomers; i++) {
		monomers[i].r.vector[0] = spacing*(double)i;
		monomers[i].r.vector[1] = 0;
		monomers[i].r.vector[2] = 0;
	}
}

int branchedChain::dendrimerMass(const int f, const int g) {
	return (2-f*round(pow( (double)(f-1), (double)g )) )/(2-f);
}

int branchedChain::arcLength(const int start, const int nofollow){
	int countedMonomers = 0;
	if(start >= numMonomers || start < 0 || nofollow >= numMonomers || nofollow < 0){
		std::cerr << "Indexing non-existant monomner"; exit(0);
	}
	//count away branches
	for(int i = 0; i < MAXBRANCHES; i++) {
		int away = monomers[start].branch[i];
		if(away != start && away != nofollow) { //away, as does not point to self
			std::cout << diffSquared( &monomers[start].r, &monomers[away].r) << ",";
			countedMonomers += this->arcLength(away,start);
		}
	}
	return countedMonomers + 1;
}

//Computes statistics about average edge length.  edgeLabels is 2x as long as edges
//This labels edged in what I call standard order.  Scan adjacency row-first.  See structure
//of double loop
void branchedChain::edgeLength(OnlineVar *el, OnlineVar *avgPos, OnlineVar *seg_lens, int *labeltrans){
	int sd;
	int en = 0;
	int sl = 0;
	jVector origin = {{0,0,0}};
	jVector t_center;
	if(numPhantoms > 0) {sd = REPEAT_PHANTOMS + 1;} else {sd = 1;}
	for(int i = 0; i < numMonomers - numPhantoms; i++){
		for(int j = i + 1; j < numMonomers - numPhantoms; j++){
			if(distMatrix[i*numMonomers + j] == sd) {
				el[en].addValue( sqrt(diffSquared(&monomers[i].r, &monomers[j].r)) ); //edge length
				vectorAvgP(&t_center, &monomers[i].r, &monomers[j].r); //find the center of the edge
				avgPos[en].addValue( sqrt(diffSquared(&t_center,&origin)) ); //find the distance from the origin 
				en++;
				int end = j;
				int middle = -666;
				int begin = -666;
				if(sd > 2) { //Find average length of phantom segments
					while( begin != i){
						middle = parentMatrix[i*numMonomers + end];
						begin = parentMatrix[i*numMonomers + middle];
						seg_lens[sl].addValue( sqrt(diffSquared(&monomers[begin].r, &monomers[end].r)) ); //edge length
						end = begin;
						labeltrans[sl] = en - 1; //I subtract 1 because en was incremented earlier.
						sl++;
					}
				}
			}
		}
	}
}

void branchedChain::insertPhantoms(const int RegularMonomers){
	//makePhantoms2(RegularMonomers,0,0,0);
	std::vector<int> q;
	makePhantoms2(RegularMonomers,0,0, &q);
	if( 2*(RegularMonomers - 1) != q.size()) {
		std::cerr << "Your graph has cycles.  Cannot handle that !" << q.size() << std::endl;
		exit(0);
	}
	for(int i = 0; i < RegularMonomers - 1; i++) {  //number of edges
		int s = i*2;  int a = s + 1;  //iterate over edges in q
		jVector diff, newPos, scaled;
		unlinkTwo(q[s],q[a]);  //unlink start and away
		linkTwo(q[s], RegularMonomers + numPhantoms); //link in the phantom

		vectorSubP(&diff, &monomers[q[a]].r, &monomers[q[s]].r );
		scaleP(&scaled, &diff, 1.0/(REPEAT_PHANTOMS+1.0) );
		vectorAddP(&monomers[RegularMonomers + numPhantoms].r, &scaled, &monomers[q[s]].r);

		int firstPhant = RegularMonomers + numPhantoms;
		for(int ii = 0; ii < REPEAT_PHANTOMS -1; ii++ ){ //link in middle phantsom
			linkTwo(RegularMonomers + numPhantoms, RegularMonomers + numPhantoms + 1);
			monomers[RegularMonomers + numPhantoms + 1].zeroBranch = RegularMonomers + numPhantoms ;  //set the correct zero branch
			vectorAddP(&monomers[RegularMonomers + numPhantoms + 1].r, &scaled, &monomers[RegularMonomers + numPhantoms].r);
			numPhantoms++;
		}
		
		linkTwo(RegularMonomers + numPhantoms, q[a]); //link in the phantom

		monomers[q[a]].zeroBranch = RegularMonomers + numPhantoms;  //set the correct zero branch
		monomers[firstPhant].zeroBranch = q[s];
		

		numPhantoms++;
	}
	buildParentMatrix();
}

int branchedChain::makePhantoms2(const int numRegularMon, const int start, const int nofollow, std::vector<int> *q){
	int countedMonomers = 0;
	if(start >= numRegularMon || start < 0 || nofollow >= numRegularMon || nofollow < 0){
		std::cerr << "Indexing non-existant monomner"; exit(0);
	}
	//count away branches
	for(int i = 0; i < MAXBRANCHES; i++) {
		int away = monomers[start].branch[i];
		if(away != start && away != nofollow) { //away, as does not point to self
			//std::cout << diffSquared( &monomers[start].r, &monomers[away].r) << ",";
			q->push_back(start); q->push_back(away);  //list of edges
			countedMonomers += this->makePhantoms2(numRegularMon, away,start, q);
		}
	}
	return countedMonomers + 1;
}

//compact way to set monomer position
void branchedChain::smPos(const int i, const double x, const double y, const double z) {
	monomers[i].r.vector[0] = x;
	monomers[i].r.vector[1] = y;
	monomers[i].r.vector[2] = z;
}

//recursively add to dendrimer
void branchedChain::_addToDend(int *cm, const int mon, const int f, const int g, const double x, const double y, const double z, stateGen *generator) {
	//rng maps
	std::uniform_real_distribution<double> phim(0,M_PI*2.0);
	std::uniform_real_distribution<double> ca(-1.0,1.0);

	if(g == 1){return;} //base condition

	for(int i=0; i < f -1; i++){
		double theta = acos(ca(*generator)); double phi = phim(*generator);
		double dx,dy,dz;
		if(space_dim == 2) { //restrict to x-y plane
			dx = spacing*cos(phi);
			dy = spacing*sin(phi);
			dz = 0;
		} else {
			dx = spacing*sin(theta)*cos(phi);
			dy = spacing*sin(theta)*sin(phi);
			dz = spacing*cos(theta);
		}
		double nx = x+dx; double ny = y+dy; double nz = z+dz;
		linkTwo(mon,*cm);
		monomers[*cm].zeroBranch = mon;
		smPos(*cm,nx,ny,nz); *cm = *cm + 1;
		_addToDend(cm,*cm-1,f,g-1,nx,ny,nz,generator);
	}

}

void branchedChain::createDendrimerR(const int mass, const int f, const int g, stateGen *generator){
	if(mass!= dendrimerMass(f,g) ) {
		std::cout << "Error: numMonomers" << mass <<" is not correct for a dendrimer:" << dendrimerMass(f,g) << std::endl;
		exit(0);
	}
	//rng maps
	std::uniform_real_distribution<double> phim(0,M_PI*2.0);
	std::uniform_real_distribution<double> ca(-1.0,1.0);
	//reset chain
	//Initialize the array
	for(int i = 0; i < mass; i++) {
		smPos(i,0,0,0);
		//these branches index an array element.
		//Link each monomer to itself by default.  gas
		for(int j = 0; j < MAXBRANCHES; j++){
			monomers[i].branch[j]=i;
		}
	}

	//initialize dendrimer's first generation and core
	int cm = 0;
	monomers[0].zeroBranch = 0;
	smPos(0,0,0,0); cm++;
	for(int i = 0; i < f; i++) {
		double theta = acos(ca(*generator)); double phi = phim(*generator);
		double x,y,z;
		if(space_dim == 3) {
			x = spacing*sin(theta)*cos(phi); y = spacing*sin(theta)*sin(phi); z = spacing*cos(theta);
		} else {	//Restrict to x-y plane
			x = spacing*cos(phi); y = spacing*sin(phi); z = 0;
		}
		smPos(cm,x,y,z);
		monomers[cm].zeroBranch = 0;
		linkTwo(0,cm); //link core to first generation
		cm++;
	}
	//Do other generations
	for(int i = 0; i < f; i++) {
		double x = monomers[i+1].r.vector[0];double y = monomers[i+1].r.vector[1];double z = monomers[i+1].r.vector[2];
		_addToDend(&cm, i+1, f, g, x,y,z, generator);
	}
	this->buildParentMatrix();
}

void _setRandPos(branchedChain *mol, unsigned int last, stateGen * generator);
void _setRandPos(branchedChain *mol, unsigned int last, stateGen * generator) {
	std::uniform_real_distribution<double> phim(0,M_PI*2.0);
	std::uniform_real_distribution<double> ca(-1.0,1.0);

	jVector res, np;

	for(int i = 0; i < MAXBRANCHES; i++) {
		int cur = mol->monomers[last].branch[i];
		if(cur != last && cur != mol->monomers[last].zeroBranch) { //assign adjacent monomers positions
			double theta = acos(ca(*generator)); double phi = phim(*generator); //create random vector
			double x = mol->spacing*sin(theta)*cos(phi); double y = mol->spacing*sin(theta)*sin(phi); double z = mol->spacing*cos(theta);
			np.vector[0] = x; np.vector[1] = y; np.vector[2] = z; np.vector[3] = 0;
			vectorAddP(&res, &np, &mol->monomers[last].r); //add new vector to prev
			mol->smPos(cur, res.vector[0], res.vector[1], res.vector[2]); //set position
			_setRandPos(mol, cur, generator); //make recursive call
		}
	}

}

//Load a .csv adjacency matrix and build it.
void branchedChain::loadAdjacency(int mass, const char *fname, stateGen * generator) {
	int32_t dtype;
	int32_t dim= 666;
	int64_t *adj;
	int count = 0;

	FILE *fp = fopen(fname, "rb");
	if(fp == NULL) {
		std::cerr << "Could not find file "<< fname << std::endl;
		exit(0);
	}
	count += fread(&dtype, sizeof(dtype), (size_t)1, fp);
	count += fread(&dim, sizeof(dim), (size_t)1, fp);
	adj = new int64_t[dim*dim];
	count += fread(adj, sizeof(adj[0]), (size_t)dim*dim, fp);
	fclose(fp);
	if (dim != mass) {
		std::cerr << "adjacency matrix is not of the correct size, " << mass << "Size is " << dim << std::endl;
		exit(0);
	}
	std::cerr << "Type is " << dtype << std::endl;
	std::cerr << "Dimension is " << dim << std::endl;
	//determine largest functionality
	int maxfunc=0;
	int tmf;
	for(int i = 0; i < dim; i++) {
		tmf = 0;
		for(int j = 0; j < dim; j++){
			if( adj[i*dim + j] == 1){tmf++;}
		}
		if(tmf > maxfunc) { maxfunc = tmf;}
	}
	if( maxfunc > MAXBRANCHES ){
		std::cerr << "Must increase MAXBRANCHES to at least " << maxfunc << std::endl;
		exit(0);
	}
	//Connect all the nodes
	int bi;
	for(int i = 0; i < dim; i++) {
		bi = 0;
		for(int j = 0; j < dim; j++){
			if( adj[i*dim + j] == 1){
				monomers[i].branch[bi] = j;
				bi++;
			}
		}
	}
	delete[] adj;
	buildParentMatrix();
	//set zero branches
	for(int i = 0; i < dim; i++) {
		monomers[i].zeroBranch = parentMatrix[i];
	}
	//set random position
	smPos(0,0,0,0); //set monomer 0 to origin
	_setRandPos(this, 0, generator);
}

branchedChain::~branchedChain() {
	//delete[] monomers;
	delete[] _fmonom16;
	delete[] parentMatrix;
	delete[] divStack;
	delete[] distMatrix;
	delete[] selfEnergyPerParticle;
	if(de3DArray != NULL) {delete[] de3DArray;}
}

branchedChain::branchedChain(const int num, const int dim){
	if(num < 2) { std::cerr << "Must have more than one monomer per branchedChain" << std::endl;
			exit(EXIT_FAILURE); }
	try {
		_fmonom16 = new jMonomer[2*num+16];
		monomers = ((uint64_t)_fmonom16 % 16) + _fmonom16;
		//monomers = new jMonomer[2*num]; //allocate twice the amount, to create the save buffer space
		selfEnergyPerParticle = new double[num];
		parentMatrix = new int[num*num]; //allocate partent matrix
		distMatrix = new int[num*num];
		saveIndex = num;
	} catch (std::bad_alloc& ba) {
		std::cerr << "Could not allocate a branched chain"
			  << ba.what() << "\n";
		exit(EXIT_FAILURE);
	}
	//Now that that is done, on with the good stuff
	space_dim = dim;
	numMonomers = num;
	numPhantoms = 0;
	energy = 0.0;
	selfEnergy = 0.0;
	spacing = 1.0;
	epsilon = 0.0;
	beta = 0;
	sigma = 0.9; //sqrt(2) is another good default
	oldSig = -1;
	oldSpacing = -1;
	acceptedMoves = 0; movesLastTimer = 0; speedFlag = 0; indicatorCount = 0; attemptedMoves = 0;

	findDerivativeDensity = 0;
	de3Ddim = 0; 
	de3DArray = NULL; //allocate space
	//Initialize the array
	for(int i = 0; i < numMonomers; i++) {
		monomers[i].r.vector[0] = 0.0;
		monomers[i].r.vector[1] = 0.0;
		monomers[i].r.vector[2] = 0.0;
		//these branches index an array element.
		//Link each monomer to itself by default.  gas
		for(int j = 0; j < MAXBRANCHES; j++){
			monomers[i].branch[j]=i;
		}
	}
	//init matrices
	for(int i = 0; i < numMonomers*numMonomers; i++){
		parentMatrix[i] = 0;
		distMatrix[i] = -1;
	}
	//allocate stack.  Should probably allocate all the programs data in one big chunk
	//do later.
	divStack = new int[numMonomers];
}

int branchedChain::countMonomers(const int start, const int nofollow){
	int countedMonomers = 0;
	if(start >= numMonomers || start < 0 || nofollow >= numMonomers || nofollow < 0){
		std::cerr << "Indexing non-existant monomner"; exit(0);
	}
	//count away branches
	for(int i = 0; i < MAXBRANCHES; i++) {
		int away = monomers[start].branch[i];
		if(away != start && away != nofollow) { //away, as does not point to self
			countedMonomers += this->countMonomers(away,start);
		}
	}
	if(this->isStackEnabled()) {
		this->pushStack(start);  //cram us on the stack

	}
	return countedMonomers + 1;
}


//return which branch has the largest monomer count
int branchedChain::largestSection(const int start) {
	int maxSection = 0;
	int theSection = 0;
	for(int i = 0; i < MAXBRANCHES; i++){
		int away = monomers[start].branch[i];
		if(away != start) { //only follow branches pointing away
			int secSize = this->countMonomers(away,start);
			if(secSize > maxSection ) { maxSection = secSize; theSection = away;}
		}
	}
	//Now we build the stack.  Put the largest section on first, including start
	this->enableStack();
	this->pushStack(start);
	this->countMonomers(theSection, start);
	for(int i = 0; i < MAXBRANCHES; i++) {
		int away = monomers[start].branch[i];
		if(away != start && away != theSection) { //only follow smaller branched
			this->countMonomers(away,start);  //Should I push on to the stack here?
		}
	}
	divStackLargest = maxSection + 1;
	this->disableStack();
	return theSection;
}


//return which branch does not contain monomer zero
int branchedChain::nonOriginSection(const int start) {
	int secSize = 0;

	this->enableStack();
	this->pushStack(start);

	for(int i = 0; i < MAXBRANCHES; i++){
		int away = monomers[start].branch[i];
		if(away != start && away != monomers[i].zeroBranch) { //only follow branches pointing away from origin
			secSize = this->countMonomers(away,start);
		}
	}
	this->disableStack();
	return secSize;
}

bool branchedChain::rotateBranch(const int start, const int branch, const int r1, const int r2, const int r3) {
	jVector base; memcpy(&base, &monomers[start].r.vector, sizeof(jVector));	//base vector
	jMatrix rot = eulerAngleMatrix(r1,r2,r3);									//rotation
	//Now rotate all the monomers in place
	this->rotateMonomer(branch,start,&rot,&base);
	return true;
}

bool branchedChain::rotateCMBranch(const int start, const int r1, const int r2, const int r3) {
	jVector base; 
	jMatrix rot = eulerAngleMatrix(r1,r2,r3);									//rotation
	//Now rotate all the monomers in place
	if(start != 0) {
		memcpy(&base, &monomers[monomers[start].zeroBranch].r.vector, sizeof(jVector));	//base vector, back up a step
		//std::cout << "Drew monomer " << start << " Rotating monomers:" << std::endl;
		this->rotateMonomer(start,monomers[start].zeroBranch,&rot,&base);
		//std::cout << std::endl;

	} else {
		//rotate all monomers
		//memcpy(&base, &monomers[start].r.vector, sizeof(jVector));	//base vector, back up a step
		//this->rotateMonomer(start,start,&rot,&base);
	}
	return true;
}


//Yeah, it is duplicate code.  But I do not want to waste time implementing a generalized
//Recursion scheme.  Maybe when I learn Lambdas in C++ I will revisit this
//returns number of monomers rotated.
int branchedChain::rotateMonomer(const int start, const int nofollow, const jMatrix *rot, const jVector *base) {
	int countedMonomers = 0;
	if(start >= numMonomers || start < 0 || nofollow >= numMonomers || nofollow < 0){
		std::cerr << "Indexing non-existant monomner"; exit(0);
	}
	//count away branches
	for(int i = 0; i < MAXBRANCHES; i++) {
		int away = monomers[start].branch[i];
		if(away != start && away != nofollow) { //away, as does not point to self
			countedMonomers += this->rotateMonomer(away,start,rot,base);
		}
	}
	//std::cout << start << ",";
	jVector rotated = {{0,0,0}};
	jVector translated = {{0,0,0}};
	vectorSubP(&translated, &monomers[start].r, base);	//translated = start-base
	transformVectorP(&rotated, rot, &translated);				//rotated=M*translated
	vectorAddP(&monomers[start].r, &rotated, base);		//add base back and store result
	return countedMonomers + 1;  //The plus 1 is for counting self
}

//Breadth first search with path reconstruction
//This is not performance critical code
void branchedChain::rotateMonomersNR(const int i, const int nofollow, const jMatrix *rot, const jVector *base) {
	std::queue<int> qu;
	int labels[CHAINLENGTH];
	for(int i = 0; i < this->numMonomers; i++) { //initialize
		labels[i] = 0;
	}
	qu.push(i);
	labels[i] = 1; //mark start as discovered

	jVector rotated = {{0,0,0}};
	jVector translated = {{0,0,0}};

	while( !qu.empty() ){
		int cur = qu.front(); qu.pop();  //deque
		vectorSubP(&translated, &monomers[cur].r, base);	//translated = start-base
		transformVectorP(&rotated, rot, &translated);
		vectorAddP(&monomers[cur].r, &rotated, base);		//add base back and store result
		for(int j = 0;  j < MAXBRANCHES; j++){
			int db = this->monomers[cur].branch[j];
			if(db != cur && labels[db] == 0 && db != nofollow) { //if not pointing at itself, and has not been searched
				qu.push(db); labels[db] = 1;
			}
		}
	}
}

void branchedChain::save(void) {
	memcpy(&monomers[numMonomers],&monomers[0],sizeof(jMonomer)*numMonomers);
}

void branchedChain::restore(void) {
	memcpy(&monomers[0],&monomers[numMonomers],sizeof(jMonomer)*numMonomers);
}

bool branchedChain::isStackEnabled(){
	if(divStackIndex >= 0 && divStackIndex < numMonomers) {return true;} else {return false;}
}
void branchedChain::enableStack(){ //When enabled, countMonomers pushes each monomer it find on the stack
	divStackIndex = 0;
}
void branchedChain::disableStack(){ //countMonomers no longer pushes on stack
	divStackIndex = -1;
}
void branchedChain::pushStack(const int val){
	if(this->isStackEnabled()){
		divStack[divStackIndex] = val;
		divStackIndex++;
	}
}

//Run MC pivot for a certrain number of steps
int branchedChain::runMC(const int steps, stateGen *generator) {
	std::uniform_real_distribution<double> azDistribution(0,M_PI*2.0); //azimuthal angle
	std::uniform_real_distribution<double> thetaDistribution(-1.0,1.0); //cos theta
	std::uniform_int_distribution<int> 	intDistribution(0, numMonomers - 1);
    std::uniform_real_distribution<double> accdist(0.0,1.0); //metropolis
    std::uniform_int_distribution<int> latticeRot(0,3);

	double newenergy;
	this->energy = this->totalLJ()+ EXTERNAL_ENERGY_COEFF*this->externalEnergy(MY_TENSION);

	for(unsigned int i = 0; i < steps; i++) {
			//generate rotation angles
			double rz1,rx1,rz2;

			//check the speed flag generated by the interrupt
			//This is for displaying periodic progresss and
			//performance metric
			if(speedFlag){
				std::ostringstream str1,str2;
				double amp = (double)(acceptedMoves - movesLastTimer)/INDICATOR_PERIOD;
				movesLastTimer = acceptedMoves;
				double totalMoves = BATCHES*LINES*DECORR_TIME;
				str1 << "AMP/Second " << amp << ". Completed "<< 100.0*(double)acceptedMoves/totalMoves << "%.  ETA " << ((totalMoves-acceptedMoves)/amp)/3600.0 << "h.  accepted/failed moves "<<acceptedMoves<<"/" << attemptedMoves << "Total Moves: " << totalMoves << std::endl;
				str2 << "Accepted Moves Per Second " << amp << ". Completed " << acceptedMoves <<" Moves." << std::endl;
				if(indicatorCount % INDICATOR_FILE_PERIOD == 0) { //how often to write progress to file
					std::ostringstream pfn;
					pid_t myPid = getpid();
					pfn << "progress_" << myPid << ".txt";
					indicatorCount = 0;
					std::ofstream lf(pfn.str(),std::ofstream::out);
					#ifdef MODE_RUN
					lf << str1.str();
					#else
					lf << str2.str();
					#endif
					lf.close();
				}
				#ifdef MODE_RUN
				std::cerr << str1.str();
				#else
				std::cerr << str2.str();
				#endif

				alarm(INDICATOR_PERIOD);
				speedFlag = 0;
				indicatorCount++;
			}

			rx1 = acos(thetaDistribution(*generator));
			rz1 = azDistribution(*generator);
			rz2 = azDistribution(*generator);
			int randMonomer = intDistribution(*generator);
			//int bigSection = mol1.largestSection(randMonomer);
			this->save();
			if(space_dim == 2){
				this->rotateCMBranch(randMonomer,rz1,0,0);
			} else {
				this->rotateCMBranch(randMonomer,rz1,rx1,rz2);
			}

			//externalEnergy is debug
	        newenergy = this->totalLJ() + EXTERNAL_ENERGY_COEFF*this->externalEnergy(MY_TENSION); //calculate initial energy of configuration

	        //Metroplois condition
	        double pAcc = std::min(1.0, exp(-(newenergy-this->energy) ) );
	        if(accdist(*generator) > pAcc) {
	        	attemptedMoves++;
	            this->restore();
	            i += -1;
				continue;
	        } else { //accept the new move
	            this->energy = newenergy;
	            acceptedMoves++;
	        }

		}
	return 0;
}

double sincPot(double r, const double sig, const double eps);
inline double sincPot(double r, const double sig, const double eps){
	double rr = r/sig;
	if(rr < 0.01) {
		return eps;
	}
	else if( rr >= 1.0) {
		return 0.0;
	}
	else{
		return eps*sin(M_PI/sig*r)/(M_PI/sig*r);
	}
}


double branchedChain::totalLJ() {
	//return false;  //disable SARW
    double epen = 0.0;  //energy penalty
	int nm = this->numMonomers - this->numPhantoms;
	if(spacing == sigma) {
		std::cerr << "Spacing cannot exactly equal sigma!  This will bug out." << std::endl;
		exit(0);
	}
	if(oldSpacing < 0) {  //check if we need to recalcuate self energy
		this->setSelfEnergy(spacing, sigma);
		oldSpacing = spacing; oldSig = sigma;
	}
	if(spacing != oldSpacing || sigma != oldSig){
		this->setSelfEnergy(spacing,sigma);
		oldSpacing = spacing; oldSig = sigma;
	}
	double d2 = 0.0;
	int al = 1;
	double tup, tup2;
	if( numPhantoms > 0) {al = REPEAT_PHANTOMS + 1;}
#ifdef FREE_ENERGY_DENSITY
	//reset interaction energy
	for(int i = 0; i < nm; i++) {
		monomers[i].lastInteraction = 0;
	}
#endif

	for(int i = 0; i < nm; i++) {
		for(int j=i+1; j < nm; j++){
			if(distMatrix[i*numMonomers + j] > al) { //ignore nearest neighbors
				d2 = diffSquared(&monomers[i].r, &monomers[j].r);
				tup = unitPot(d2,sigma);
				epen += tup;  //These 3 lines are for energy density.  They cost a little, 10% speed
#ifdef FREE_ENERGY_DENSITY
				monomers[i].lastInteraction += tup;  //make sure to divide by 2 in post processing.  We double count
				monomers[j].lastInteraction += tup;
#endif
			}
		}
	}

	//Now find total bending energy
	jVector l1, l2;
	int connectingMonomer;
	double bendingEnergy = 0.0;
	double ls = spacing*spacing/( (REPEAT_PHANTOMS+1.0)*(REPEAT_PHANTOMS+1.0) );
	if(numPhantoms > 0){
		for(int i = 0; i < nm; i++) {
			for(int j=i+1; j < nm; j++){
				if(distMatrix[i*numMonomers + j] == REPEAT_PHANTOMS + 1) { //adjacnt pair
					int termMon = j;
					int startMon = i;
					for(int iy = 0; iy < (REPEAT_PHANTOMS+1)/2; iy ++) { //find bending energy bet
						connectingMonomer = parentMatrix[i*numMonomers + termMon]; //First connecting monomer
						startMon = parentMatrix[i*numMonomers + connectingMonomer]; //start monomer of segment
						vectorSubP(&l1, &monomers[connectingMonomer].r, &monomers[startMon].r);
						vectorSubP(&l2, &monomers[termMon].r, &monomers[connectingMonomer].r);
						bendingEnergy += std::abs(dotP(&l1, &l2)/ls);  //1.0 - 0.5*(1.0 + dotP(&l1, &l2)/ls);
						termMon = startMon;
					}
					/*connectingMonomer = parentMatrix[i*numMonomers + j]; //the monomer between i and j
					vectorSubP(&l1, &monomers[connectingMonomer].r, &monomers[i].r);
					vectorSubP(&l2, &monomers[j].r, &monomers[connectingMonomer].r);
					bendingEnergy += std::abs(dotP(&l1, &l2)/ls);  //1.0 - 0.5*(1.0 + dotP(&l1, &l2)/ls);*/
				}
			}
		}
	}
	return this->epsilon*(epen-this->selfEnergy)+beta*bendingEnergy;
}

double branchedChain::setSelfEnergy(const double spacing, const double sig) {
	//return false;  //disable SARW
    double epen = 0.0;  //energy 
	int nm = this->numMonomers - numPhantoms;
	
	buildParentMatrix();

    std::cout << "Setting self energy" << std::endl;


	int al = 1; //graph distance that counts as nearest neighbor
	double tr, acc;
	if(numPhantoms > 0){ al= REPEAT_PHANTOMS + 1; }

	for(int i = 0; i < nm; i++){
		selfEnergyPerParticle[i] = 0;
	}

	if(epsilon == 0.0) { return 0.0; }
	for(int i = 0; i < nm; i++) {
		acc = 0; //accumumlated self energy for particle i
		for(int j=i+1; j < nm; j++){
			if(distMatrix[i*numMonomers + j] > al ){ //ignore nearest neighbor interactions
				double d2 = distMatrix[i*numMonomers + j]*distMatrix[i*numMonomers + j];  //3 lines:  Find the max distance between 2 monomers
				if(numPhantoms > 0) { d2 = d2*spacing*spacing/( (REPEAT_PHANTOMS+1)*(REPEAT_PHANTOMS + 1) ); }// distances over over estimated with phantom monomers
				else { d2 = d2*spacing*spacing; }
				tr = unitPot(d2,sig);
				epen += tr;//sincPot(sqrt(d2), sig, eps);
				acc += tr;
			}
		}
		selfEnergyPerParticle[i] = acc; //set self energy per particle in array
	}
	std::cerr << "Self Energy is " << epen << std::endl;
	this->selfEnergy = epen;
	return epen;
}

//fast
double branchedChain::totalLJFast() {
	double epen = 0.0;  //energy penalty
    double dr2[sizeof(jVector)/sizeof(double)+16];
    double *d2p = (double *)__builtin_assume_aligned( (double *)( (uint64_t)dr2 % 16 + (uint64_t)dr2 ), 16);
	int nm = this->numMonomers - this->numPhantoms;
	if(spacing == sigma) {
		std::cerr << "Spacing cannot exactly equal sigma!  This will bug out." << std::endl;
		exit(0);
	}
	if(oldSpacing < 0) {  //check if we need to recalcuate self energy
		this->setSelfEnergy(spacing, sigma);
		oldSpacing = spacing; oldSig = sigma;
	}
	if(spacing != oldSpacing || sigma != oldSig){
		this->setSelfEnergy(spacing,sigma);
		oldSpacing = spacing; oldSig = sigma;
	}
	//jMonomer *monomers = this->monomers;
	for(int i = 0; i < nm; i++) {
		for(int j=i+1; j < nm; j++){
			double d2 = diffSquared(&monomers[i].r, &monomers[j].r);
			double di = sqrt(d2);
			int ns = (int)floor( (di-sigma)/spacing );
			if(ns > 1) {
				j += ns;
				continue;
			}
			epen += unitPot(d2,sigma);//sincPot(sqrt(d2), sig, eps);
		}
	}
	return this->epsilon*(epen-this->selfEnergy);
}

//Breadth first search with path reconstruction
//This is not performance critical code
void branchedChain::bfs(const int i) {
	std::queue<int> qu;
	int dist = 0;
	int *labels = new int[numMonomers];
	for(int i = 0; i < this->numMonomers; i++) { //initialize
		labels[i] = 0;
	}
	qu.push(i);
	labels[i] = 1; //mark start as discovered
	distMatrix[i*numMonomers + i] = 0; //i is 0 distance from itself
	while( !qu.empty() ){
		int cur = qu.front(); qu.pop();  //deque
		for(int j = 0;  j < MAXBRANCHES; j++){
			int db = this->monomers[cur].branch[j];
			if(db != cur && labels[db] == 0) { //if not pointing at itself, and has not been searched
				qu.push(db); labels[db] = 1;
				distMatrix[i*numMonomers + db] = distMatrix[i*numMonomers + cur] + 1;
				this->parentMatrix[i*this->numMonomers + db] = cur;
			}
		}
	}
	delete[] labels;
}

void branchedChain::buildParentMatrix() {
	for(int i = 0; i < this->numMonomers; i++){
		this->bfs(i);
	}
}

typedef struct myQ {
	int db;
	int hops;
} Skip;

/*Simple quadratic
inline double unitPot(const double d2, const double sig){
	double x = d2/(sig*sig);
	if(x < 1.0) {
		return 1.0 - x;
	} else {
		return 0.0;
	}

}*/

//fast exponential approximation
inline
double exp1(double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}

inline double unitPot(const double d2, const double sig){
	double x = d2/(sig*sig);
	return (1.0/(M_PI*sig*sig))*exp1(-x);
}

//enable the finding of derivative density.  Dim is the number of divisions per dimension
//domain is the distance from the origin the box should extend along the axis.
//This takes memory on order dim^3, so be careful
//Also call this AFTER phantoms are inserted
void branchedChain::enableDerivativeDensity(const double dim, const double domain) {
	if(findDerivativeDensity == 1) { return; }
	findDerivativeDensity = 1;
	de3Ddim = dim; 
	de3Ddomain = domain;
	unsigned int id = (unsigned int)round(dim);
	de3DArray = new double[id*id*id]; //allocate space
	for(int i = 0; i < id*id*id; i++) {
		de3DArray[i] = 0.0;
	}
}

//Takes the epsilon value for the derivative.  Returns the total derivative.  This should equal the 
//system derivative
double branchedChain::binParticles(const double de) {
	unsigned int k[3];
	double r[3];
	double *x;
	double binWidth = 2.0/de3Ddim;
	unsigned int id = (int)round(de3Ddim);
	int *stack;
	double totalDE = 0;

#ifndef FREE_ENERGY_DENSITY
	std::cerr << "You must enable the FREE_ENERGY_DENSITY directive to use this function! " << std::endl;
	exit(0);
#endif

	if(findDerivativeDensity != 1) {
		std::cerr << "Did not initialize energy density!" << std::endl;
		exit(0);
	}

	//reset 3D grid
	for(int i = 0; i < id*id*id; i++) {
		de3DArray[i] = 0.0;
	}

	for(int j = 0; j < numMonomers - numPhantoms; j++){
		x = monomers[j].r.vector;
		for(int i = 0; i < 3; i++) { r[i] = 1.0 + x[i]/ de3Ddim;} //shift and rescale x coords to r
		for(int i = 0; i < 3; i++) { if(r[i] < 0.0 || r[i] >= 2.0 ) { continue;} } //check rescaling
		for(int i = 0; i < 3; i++) { k[i] = (unsigned int)floor(r[i]/binWidth);  }  //convert to 3D index
		double baseEnergy = monomers[j].lastInteraction/2.0 - selfEnergyPerParticle[j]; //the divide by 2 is we overcount in totalLJ
		double dec = ((this->epsilon + de)*baseEnergy - this->epsilon*baseEnergy)/de; //freeEnergy derivative for this particle
		totalDE +=  dec;
		unsigned int indx = k[0]*id + k[1] + id*id*k[2];
		if( indx >= id*id*id) {
			std::cerr <<" Array out of bounds!" << de3Ddim << std::endl;
			for(int p = 0; p < 3; p++){
				std::cerr << "k["<<p<<"] is " << k[p] << std::endl;
			}
		}
		de3DArray[k[0]*id + k[1] + id*id*k[2]] += dec; //add particle to cell in 3D
	}
	return totalDE;
}

double branchedChain::energyDerivative(const double de){
	double old_energy = energy;
	double old_eps = this->epsilon;
	this->epsilon += de;
	double new_energy = totalLJ();
	double derp = (new_energy - old_energy)/de;
	this->epsilon = old_eps;
	return derp;
}

double branchedChain::externalEnergy(const double T) {
	double r = sqrt(diffSquared(&monomers[0].r, &monomers[numMonomers- numPhantoms - 1].r));
	return -T*r;
}


void branchedChain::findDensity2D(Eigen::Matrix3Xd pos, double *dens, const int bins, const double rmax){
	jVector rcm, monomer;
	Eigen::Matrix3Xd rcmv = pos.rowwise().sum()/(double)pos.cols(); //finds CM
	for(int d = 0; d < 3; d++) { //Copy Eigen vector to my jVector
		rcm.vector[d] = rcmv(d,0);
	}
	for(int i = 0; i < numMonomers - numPhantoms; i++){
		jVector d = {{0,0,0}};
		for(int d = 0; d < 3; d++){monomer.vector[d] = pos(d, i); } //copy to monomer
		vectorSubV(&d,&monomer, &rcm); // Monomer offset from CM
		int pbx = (int) floor(bins/2.0*(1.0 + d.vector[0]/rmax) );	// x bin
		int pby = (int) floor(bins/2.0*(1.0 + d.vector[1]/rmax) );	// y bin
		if(pbx >= bins || pby >= bins || pbx < 0 || pby < 0) {
			std::cout << "Your rmax is too small!" <<std::endl;
			exit(0);
		}
		dens[pbx + bins*pby] += 1.0; //normalization such that integra lof density is 1
	}
}

void branchedChain::findDensity2D(double *dens, const int bins, const double rmax){
	jVector rcm;
	findCm(&rcm);
	for(int i = 0; i < numMonomers - numPhantoms; i++){
		jVector d = {{0,0,0}};
		vectorSubV(&d,&monomers[i].r, &rcm); // Monomer offset from CM
		int pbx = (int) floor(bins/2.0*(1.0 + d.vector[0]/rmax) );	// x bin
		int pby = (int) floor(bins/2.0*(1.0 + d.vector[1]/rmax) );	// y bin
		if(pbx >= bins || pby >= bins || pbx < 0 || pby < 0) {
			std::cout << "Your rmax is too small!" <<std::endl;
			exit(0);
		}
		dens[pbx + bins*pby] += 1.0; //normalization such that integra lof density is 1
	}
}

void branchedChain::findDensity(double *dens, const int bins, const double rmax){
	jVector rcm;
	findCm(&rcm);
	for(int i = 0; i < bins; i++) {	dens[i] = 0; }
	for(int i = 0; i < numMonomers - numPhantoms; i++){
		double d = sqrt(diffSquared(&monomers[i].r, &rcm) );
		int pb = (int) floor(bins*d/rmax);
		if(pb > bins) {
			std::cout << "Your rmax is too small!" <<std::endl;
			exit(0);
		} 
		dens[pb] += 1.0; //normalization such that integra lof density is 1
	}
}

void branchedChain::findDensitySingle(int mon, jVector *rcm, double *dens, const int bins, const double rmax){
	for(int i = 0; i < bins; i++) {	dens[i] = 0; }
	double d = sqrt(diffSquared(&monomers[mon].r, rcm) );
	int pb = (int) floor(bins*d/rmax);
	if(pb > bins) {
		std::cout << "Your rmax is too small!" <<std::endl;
		exit(0);
	}
	dens[pb] += 1.0; //normalization such that integra lof density is 1
}

//check if two monomers are topologically neirest neighbohrs
int branchedChain::isNearestNeighbor(const int i, const int j) {
	if(i == j) {  //also include self a NN
		return 1;
	}
	for(int k = 0; k < MAXBRANCHES; k++) {
		if(this->monomers[i].branch[k] == j || this->monomers[j].branch[k] == i) {
			return 1;
		}
	}
	return 0;
}

//write state information to file
int branchedChain::writeHeader(dataFP fp, stateGen *gen) {
	int count = 0;
	uint32_t head = 0xDEADBEEF; //DEDBEEF
	count += fwrite((void *)&head, sizeof(head),(size_t)1,fp);
	count += fwrite((void *)&gen->Seed, sizeof(gen->Seed), (size_t)1, fp); //write Seed
	count += fwrite((void *)&gen->zCount, sizeof(gen->zCount), (size_t)1, fp); //write the number of calls to discard
	count += fwrite((void *)&numMonomers, sizeof(numMonomers), (size_t)1, fp); //write chain length
	return count;
}

//write parameters to file
int branchedChain::writeParameterLine(dataFP fp) {
	uint8_t head = (uint8_t)'P';
	int count = 0;
	count += fwrite((void *)&head, sizeof(head), (size_t)1, fp);
	count += fwrite((void *)&spacing, sizeof(spacing), (size_t)1, fp); //write spacing
	count += fwrite((void *)&sigma, sizeof(sigma), (size_t)1, fp);
	count += fwrite((void *)&epsilon, sizeof(epsilon), (size_t)1, fp); //write chain length
	return count;
}

int branchedChain::writeDataLine(dataFP fp, stateGen *gen) {
	uint8_t head = (uint8_t)'D';
	int count = 0;
	count += fwrite((void *)&head, sizeof(head), (size_t)1, fp); //write line ID
	count += fwrite((void *)&gen->zCount, sizeof(gen->zCount), (size_t)1, fp); //write zCount
	//write out all monomers x,y,z position
	for(int i = 0; i < numMonomers; i++) {
		count += fwrite((void *)&monomers[i].r, sizeof(monomers[0].r), (size_t)1, fp);
	}
	return count;
}

//read state information from file
int branchedChain::readLine(dataFP fp, stateGen *gen) {
	uint8_t head;
	int count = 0;
	stateGen ts;
	count += fread((void *)&head, sizeof(head), (size_t)1, fp); //read line type, "P" or "D"
	if((char)head == 'P') {
		count += fread((void *)&spacing, sizeof(spacing), (size_t)1, fp);
		count += fread((void *)&sigma, sizeof(sigma), (size_t)1, fp);
		count += fread((void *)&epsilon, sizeof(epsilon), (size_t)1, fp);
		std::cout << spacing << sigma << epsilon << std::endl;
	}else if ((char)head == 'D') {
		count += fread((void *)&ts.zCount, sizeof(ts.zCount), (size_t)1, fp); //read zCount
		gen->seed(gen->Seed); //re-initialize generator
		gen->discard(ts.zCount); gen->zCount = ts.zCount;
		//load monomers
		for(int i = 0; i < numMonomers; i++) {
			count += fread((void *)&monomers[i].r, sizeof(monomers[0].r), (size_t)1, fp);
		}
	} else {
		std::cerr << "Invalid data line!  Corrupt file format!" << std::endl;
		exit(0);
	}
	return 1;
}

//read state information from file
int branchedChain::readHeader(dataFP fp, stateGen *gen) {
	uint32_t head;
	stateGen ts;
	fread((void *)&head, sizeof(head), (size_t)1, fp); //write Seed
	if(head == 0xDEADBEEF) { 
		fread((void *)&gen->Seed, sizeof(gen->Seed), (size_t)1, fp); //write Seed
		fread((void *)&ts.zCount, sizeof(ts.zCount), (size_t)1, fp); //write the number of calls to discard
		fread((void *)&numMonomers, sizeof(numMonomers), (size_t)1, fp); //write the number of calls to discard
		gen->seed(gen->Seed); //re-initialize generator
		gen->discard(ts.zCount); gen->zCount = ts.zCount;
	} else {
		std::cerr << "Input file is not correct format" << std::endl; 
	}
	return 1;
}

void branchedChain::saveState(const char *file, stateGen *gen, int *batch){
	dataFP fp;
	fp = fopen(file,"wb");
	if(ferror(fp)) {
		std::cerr << "Could not open file " << file <<" for writing" <<std::endl;
		exit(0);
	}
	fwrite((void *)batch, sizeof(*batch), (size_t)1, fp); //write batch number
	writeHeader(fp,gen);
	writeParameterLine(fp);
	writeDataLine(fp,gen);
	fclose(fp);
}

void branchedChain::loadState(const char *file, stateGen *gen, int *batch){
	dataFP fp;
	fp = fopen(file,"rb");
	if(ferror(fp)) {
		std::cerr << "Could not open file " << file <<" for Reading" <<std::endl;
		exit(0);
	}
	fread((void *)batch, sizeof(*batch), (size_t)1, fp);
	readHeader(fp,gen);
	readLine(fp,gen);
	readLine(fp,gen);
	fclose(fp);
}

void branchedChain::printParams(std::ostream& os){
	os << "Simulation parameters" << std::endl;
	os << "Spacing: " << spacing << std::endl;
	os << "Sigma: " << sigma << std::endl;
	os << "Epsilon: " << epsilon << std::endl;
}

Eigen::Matrix3Xd branchedChain::getMonomerPositions(){
	int nm = numMonomers - numPhantoms;
	Eigen::Matrix3Xd in(3, nm);
	for (int row = 0; row < in.rows(); row++) {
		for (int col = 0; col < in.cols(); col++) {
			in(row, col) = monomers[col].r.vector[row];
		}
	}
	return in;
}

RANDOMCLASS::result_type stateGen::operator()() {
	zCount++;
	return RANDOMCLASS::operator()();
}

void stateGen::seed( RANDOMCLASS::result_type value) {
	zCount = 0;
	Seed = value;
	RANDOMCLASS::seed(value);
}

stateGen::stateGen() {
	seed(RANDOMCLASS::default_seed);
}
