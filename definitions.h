/*
 * definitions.h
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define USE_EXTERNAL_PARAMS

//#define RUN_TESTS  //Enable to run tests


#ifdef USE_EXTERNAL_PARAMS
	#include "temp_params.h"
#else
	#define INIT_SEED 5757676576
	#define ZIPNAME "Linear400.zip"
	#define BATCHES 1//0
	#define LINES 100000//00
	#define TOPOLOGY_FILE "topo.csv"
	#define MOL_EPSILON 0.0 //Molecule parameter
#endif

#ifndef SPACE_DIMENSION
	#define SPACE_DIMENSION 3
#endif

#define LOG_STATE 10

 //Molecule parameters
 #define MOL_SPACING 1.0/0.4 //2.0/1.39388
 #define MOL_SIGMA	 0.5//MOL_SPACING*0.95
 #define MOL_BETA	 4.0 //1.0  //Bending energy

//External field
#define EXTERNAL_ENERGY_COEFF 0.0
#define MY_TENSION 0.0

#define TI_FILE "logfile.csv" //Thermodynamic integration state file
#define ADJ_FILE "Adjacency.dat"

//#define FREE_ENERGY_DENSITY  //enable energy density calculations.  Slower

//#define MODE_RELAX //the run mode.  Only one of these should be uncommented at a time
#define MODE_RUN

//The maximum functionality of the monomer.  Anything less than 2 is a gas
#define DENDRIMER
//#define LINEAR
//#define ARBITRARY

#ifndef FUNCTIONALITY  //this if block is in the event FUNCTIONALITY is specified in "temp_params.h"
	#define FUNCTIONALITY 3
#endif 
#define REPEAT_PHANTOMS 5 //#Must be Odd! 0,1,3,5,7,9,...
#define PHANTOMS 9*REPEAT_PHANTOMS
#define CHAINLENGTH 10 + PHANTOMS
#define DECORR_TIME 2500//0//0//0//0 //number of moves to relax the energy
#define EQUILIBRATION DECORR_TIME//*10*5
#define SNAPSHOT DECORR_TIME/10
//#define MONOMER_RADIUS 1.1//1.414213562
//#define EPSILON 5.0 //2

#define DE 0.001 //derivative increment
#define DENSITY_BINS 1000  //changing this will require manually chaning the python analysis file
#define DENSITY_MAX_DOMAIN 5.0 //CHAINLENGTH*MOL_SPACING //This is probably too conservative

//2D density uses the bounds above
#define DENSITY_BINS_2D 128 //Use power of 2

 //Low Level.  Do not change unless you know what you are doing
 #define MAXBRANCHES FUNCTIONALITY //This is low level.  Only change FUNCTIONALITY

 //Just Run mode parameters
 #define FRAME_PREFIX "Dend5g"

 #define INDICATOR_PERIOD 5.0
 #define INDICATOR_FILE_PERIOD 10 //how often to write indicator to a file


//These are for calculating monomer histograms.  This should be disabled by default
 /*#define NUMBINS 500
 #define HISTLEFT -10.0
 #define HISTRIGHT 10.0
 #define HISTMONOMER 1  //index of monomer to find distribution of*/

#endif /* DEFINITIONS_H_ */
