/*
 * branchedChain.h
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */

#ifndef BRANCHEDCHAIN_H_
#define BRANCHEDCHAIN_H_

#include <random>
#include "data.h"
#include <stdio.h>
#include <stdlib.h>
#include "onlineVariance.h"
#include "Kabsch.h"

#define RANDOMCLASS std::mt19937_64
typedef FILE* dataFP;

class stateGen :  public RANDOMCLASS {
public:
	uint32_t zCount;
	uint64_t Seed;
	RANDOMCLASS::result_type operator()();
	void seed( RANDOMCLASS::result_type value = RANDOMCLASS::default_seed );
	stateGen();
};

class branchedChain {
	public:
		jMonomer *monomers; 			//array of monomers
		int *parentMatrix;				//array to reconstruct paths
		int *distMatrix; 				//distances
		int    numMonomers; 			//duh
		int    numPhantoms; // the number of phantom monomers
		int space_dim; //spatial dimensions.  Options are 2 and 3
		double energy; //the energy of the system
		double selfEnergy; //self energy of nearest neighbor interactions (to be subtracted off)
		double oldSpacing, oldSig;
		double spacing;
		double sigma; //monomer interaction distance
		double epsilon; //interaction energy
		double beta; //bending energy

		unsigned long int acceptedMoves, attemptedMoves, movesLastTimer;
		int speedFlag; int indicatorCount;
		/*
		 * I implement a stack.  The purpose is to store a list of which mononmers are
		 * in the small half, and which in the large half, when I select a pivot monomer
		 * The reason I need this is because when I check for interactions, I only need
		 * to check monomers on different halves, assuming hard spheres, once in
		 * equilibrium.  The pivot move obviously will not cause overlap of monomers
		 * in the same half.
		 */
		int divStackIndex;				//index in to stack.  -1 means disabled
		int divStackLargest;			//beginning of the smaller chunk (a stack index). "Division point"
		int *divStack;					//stack to store halves of pivot molecules
		//branchedChain(const int num, jMonomer *data, int *stackspace);
		branchedChain(const int num, const int dim);	//number of monomers in molecule
		~branchedChain(); //destructor
		void saveTopofile(char *fname);
		void loadAdjacency(int mass, const char *fname, stateGen * generator);
		void createLinear(int cl);			//link the monomers in a linear chain
		int dendrimerMass(const int f, const int g); //number of monomers in a dendrimer of g generations, f functionality
		void _addToDend(int *cm, const int mon, const int f, const int g, const double x, const double y, const double z, stateGen *generator);
		void createDendrimerR(const int mass, const int f, const int g, stateGen *generator); //recursive
		int arcLength(const int start, const int nofollow);
		void edgeLength(OnlineVar *el, OnlineVar *avgDist, OnlineVar *seg_lens, int *labeltrans);
		void insertPhantoms(const int RegularMonomers);
		int makePhantoms2(const int numRegularMon, const int start, const int nofollow, std::vector<int> *q);
		void printConnections();
		void setLinearPositions(const double spacing);
		void save();					//Saves the current monomer branches and positions
		void restore();					//restores the monomer branches and positions
		/*count the number of monomers, starting at start (inclusive)
		 * and do not follow paths back to nofollow.  Nofollow is almost
		 * always start
		 */
		int countMonomers(const int start, const int nofollow); //Count monomers starting at start, not following nofollow
		int largestSection(const int start);	//Returns the branch (monomer) away from start that has the most monomers
		int nonOriginSection(const int start);
		//Higher level Methods:
		bool rotateBranch(const int start, const int branch, const int r1, const int r2, const int r3); //rotates branch about start
		bool rotateCMBranch(const int start, const int r1, const int r2, const int r3); //do not rotate zero monomer
        void enableStack(); //When enabled, countMonomers pushes each monomer it finds on the stack
		void disableStack(); //countMonomers no longer pushes on stack
		void smPos(const int i, const double x, const double y, const double z); //set monomer position

		double findRg(); //find radius of gyration]]
		double findRg(const int dim);

		//return monomer positions as eigen matrix
		Eigen::Matrix3Xd getMonomerPositions();

		void findCm(jVector *res); //finds the center of mass
		void findStructureFactorRadial(double *result, size_t rlen, double qmax);
		double externalEnergy(const double coeff);
		double totalLJ();
		double setSelfEnergy(const double spacing, const double sig);
		double totalLJFast(); //linear polymer only
		int runMC(const int steps, stateGen *generator);
		double energyDerivative(const double de);

		void findDensity2D(Eigen::Matrix3Xd pos, double *dens, const int bins, const double rmax);
		void findDensity2D(double *dens, const int bins, const double rmax);

		void findDensity(double *dens, const int bins, const double rmax); //finds density (spherical coordinates)
		void findDensitySingle(int mon, jVector *rcm, double *dens, const int bins, const double rmax); //find density of individual monomer

		void enableDerivativeDensity(const double dim, const double domain);
		double binParticles(const double de); //place particles in bin

		void bfs(const int i); //breadth first search of topology
		int isNearestNeighbor(const int i, const int j);
		void buildParentMatrix();

		int writeHeader(dataFP fp, stateGen *gen); //writes binary header
		int writeParameterLine(dataFP fp);
		int writeDataLine(dataFP fp,stateGen *gen);
		int readHeader(dataFP fp, stateGen *gen);
		int readLine(dataFP fp,stateGen *gen);

		void saveState(const char *file, stateGen *gen, int *batch);
		void loadState(const char *file, stateGen *gen, int *batch);

		void printParams(std::ostream& os);

	private:
		int saveIndex;					//the index in to monomers where the save data is
		bool isStackEnabled();
		void pushStack(const int val);
		//Link two monomers, unlinking all branches
		void linkTwo(const int src, const int dest);	//Links two monomers
		void unlinkTwo(const int src, const int dest);
		int rotateMonomer(const int start, const int nofollow, const jMatrix *rot, const jVector *base); //recursive in-place monomer rotation
		void rotateMonomersNR(const int i, const int nofollow, const jMatrix *rot, const jVector *base);
		jMonomer *_fmonom16;
		double *selfEnergyPerParticle;
		//The following are for binning particles
		unsigned int findDerivativeDensity;
		double de3Ddim, de3Ddomain;
		double * de3DArray;
};



#endif /* BRANCHEDCHAIN_H_ */
