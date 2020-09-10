/*

===============================================================
HCMM mobility model
					
Copyright (C) Chiara Boldrini, chiara.boldrini@iit.cnr.it

Based on a Community Based Mobility Model, Copyright (C) Mirco 
Musolesi University College London, m.musolesi@cs.ucl.ac.uk

===============================================================


This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
USA

###
*/


#define SOCIAL_WEIGHTS 0		// 1 for printing the cell attractivity
#define FREQUENTIAL_WEIGHTS 0	// 1 for printing the visiting frequency
#define DISTANCE_DISTRIBUTION 0 // 1 for printing the distances covered by the nodes
#define CONTROLLED_REWIRING 1 	// 1 for a single intercommunity node

#ifdef WIN32
#include "getopt.h"
#include <windows.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <istream>
#include <iostream>
#include <cstdlib>
#include <ios>
#include <sstream>		
#include <vector>
#include <map>
#include "socialNet.h"
#include "xmlParser.h"	
using namespace std;

#define RAND_INT(l,h) (((double)(random() * ((double)(h)-(l)+1.0))) + (l))

#define ALLOWED_NEIGH 1

struct hostsItem {
	double currentX;
	double currentY;
	double relativeX;
	double relativeY;
	double goalRelativeX;
	double goalRelativeY;
	double goalCurrentX;
	double goalCurrentY;
	double previousGoalX;
	double previousGoalY;
	int cellIdX;
	int cellIdY;
	double speed;
	double absSpeed;
	double af;
	bool isATraveller;
};

struct cellsItem {
	double minX;
	double minY;
	int numberOfHosts;
};

bool imInMyHome (int nowX, int nowY, int homeX, int homeY);

void parsingConfiguration (int &numHosts, double &minHostSpeed,
			   double &maxHostSpeed,
			   double &connectionThreshold,
			   double &sideLength_WholeAreaX,
			   double &sideLength_WholeAreaY, int &numberOfRows,
			   int &numberOfColumns, double &radius,
			   double &totalSimulationTime,
			   double &reconfigurationInterval,
			   double &rewiringProb, int &numberOfGroups,
			   int &inputSeed, int &numberOfTravellers,
			   double &travellerSpeed, double &travProb, double &remainingProb,
			   string fileName, vector < bool > getOpt);

/* linear congruential generator.  Generator x[n+1] = a * x[n] mod m */

//default seed 
static unsigned SEED = 93186752;

double
getRandomDouble ()
{

/* input parameters*/

/* static int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901; */

	static unsigned int a = 1223106847, m = 4294967291U, q = 3, r = 625646750;

/* static unsigned int a = 279470273, m = 4294967291U, q = 15, r = 102913196;*/
/* static unsigned int a = 1583458089, m = 2147483647, q = 1, r = 564025558; */
/* static unsigned int a = 784588716, m = 2147483647, q = 2, r = 578306215;  */
/* static unsigned int a = 16807, m = 2147483647, q = 127773, r = 2836;	  */
/* static unsigned int a = 950706376, m = 2147483647, q = 2, r = 246070895;  */

	SEED = a * (SEED % q) - r * (SEED / q);
	return ((double) SEED / (double) m);
}


// Sets the random seed of the random number generator
void
setSeed (unsigned int init)
{
	if (init != 0.0)
	SEED = init;
}


// Returns a randomly generated double number in a specific range
double
generateRandomDouble (double min, double max)
{
//double number=getRandomDouble() * (max-min+1.0) + min;
	double number = getRandomDouble () * (max - min) + min;

//double number=min+((max-min)*getRandomDouble()/RAND_MAX);
	return number;
}

// Returns a randomly generated integer number in a specific range
int
generateRandomInteger (int min, int max)
{

	int number =
	(int) floor (generateRandomDouble ((double) min, (double) max) + 1);
//printf("%s %d\n","The random generated integer number is",number);

	return number;
}


// Print on stdout the header
void
header ()
{
	printf ("%s", "\nMobility Patterns Generator for ns-2 simulator\n");
	printf ("%s", "Based on the Community Based Mobility Model\n");
	printf ("%s", "Version 0.3 December 2006\n");
	printf ("%s", "Copyright (C) Mirco Musolesi University College London\n");
	printf ("%s", "  m.musolesi@cs.ucl.ac.uk\n");
}


// Print on stdinput the usage
void
usage ()
{
	printf ("%s", "Usage: movGen [-options] \n");
	printf ("%s", "Input parameters:\n");
	printf ("%s", "-h  shows help\n");
	printf ("%s", "-x  generates XML\n");
	printf ("%s", "-n  sets the number of hosts\n");
	printf ("%s", "-t <totalSimulationTime>	 sets the total simulation time\n");
	printf ("%s", "-r <reconfigurationInterval> sets the reconfiguration interval\n");
	printf ("%s", "-s <lowerBoundSpeedHost>	 sets lower bound of the speed of the hosts (in m/s)\n");
	printf ("%s", "-S <upperBoundSpeedHost>	 sets upper bound of the speed of the hosts (in m/s)\n");
	printf ("%s", "-p <connectionThreshold>	 sets the connection threshold\n");
	printf ("%s", "-X <sideLengthXcoordinates>  sets the side length of the the simulation area - x coordinates\n");
	printf ("%s", "-Y <sideLengthYcoordinates>  sets the side length of the the simulation area - y coordinates\n");
	printf ("%s", "-R <numberOfRows> sets the number of rows\n");
	printf ("%s", "-C <numberOfColumns> sets the number of columns\n");
	printf ("%s", "-T <transmissionRange>    sets the transmission range (in m)\n");
	printf ("%s", "-w <rewiring probability> sets the rewiring probability\n");
	printf ("%s", "-G <numberOfGroups>  sets the initial number of groups for the Caveman model\n");
	printf ("%s", "-g <seedRNG>  sets the seed of the Random Number Generator\n");
	printf ("%s", "-c <numberOfTravellers> sets the number of travellers\n");
	printf ("%s", "-v <travellersSpeed> sets the speed of the travellers\n");
	printf ("%s", "-a <on/off>  sets the Girvan-Newman algorithm on/off\n");
	printf ("%s", "-d <on/off>  sets the deterministic selection of the nodes on/off. If true the selection is deterministic, if false is probabilistic\n");
	printf ("%s", "-A <on/off>  sets the colocation traces on/off\n");
	printf ("%s", "-b <on/off>  sets the communities traces on/off\n");
}


// Generates XML header 
void
generateXMLHeader (FILE * outputFileXML, double sideLength_WholeAreaX,
		   double sideLength_WholeAreaY, double transmissionRange,
		   int numberOfNodes)
{

	fprintf (outputFileXML, "<?xml version=\"1.0\"?>\n");
	fprintf (outputFileXML,
		 "<simulation xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
	fprintf (outputFileXML,
		 "xsi:noNamespaceSchemaLocation=\"http://www.i-u.de/schools/hellbrueck/ansim/xml/simulation.xsd\">\n");
	fprintf (outputFileXML, "\t<parameter>\n");
	fprintf (outputFileXML, "\t\t<area_shape>rectangle</area_shape>\n");
	fprintf (outputFileXML, "\t\t<xsize>%f</xsize>\n", sideLength_WholeAreaX);
	fprintf (outputFileXML, "\t\t<ysize>%f</ysize>\n", sideLength_WholeAreaY);
	fprintf (outputFileXML, "\t\t<extension_factor>1</extension_factor>\n");
	fprintf (outputFileXML, "\t\t<numberOfNodes>%d</numberOfNodes>\n",
		 numberOfNodes);
	fprintf (outputFileXML, "\t\t<range>%f</range>\n", transmissionRange);
	fprintf (outputFileXML,
		 "\t\t<mobility_model>socialFoundedMM</mobility_model>\n");
	fprintf (outputFileXML, "\t</parameter>\n");

}


// Generates XML section related to the definition of the node position
void
generateXMLSetNodePosition (FILE * outputFileXML, int node, double x,
				double y, double sideLength_WholeAreaX,
				double sideLength_WholeAreaY)
{
	fprintf (outputFileXML, "\t\t<node>\n");
	fprintf (outputFileXML, "\t\t\t<node_id>%d</node_id>\n", node);
	fprintf (outputFileXML, "\t\t\t<position>\n");
	fprintf (outputFileXML, "\t\t\t\t<xpos>%f</xpos>\n",
		 x - sideLength_WholeAreaX / 2.0);
	fprintf (outputFileXML, "\t\t\t\t<ypos>%f</ypos>\n",
		 y - sideLength_WholeAreaY / 2.0);
	fprintf (outputFileXML, "\t\t\t</position>\n");
	fprintf (outputFileXML, "\t\t</node>\n");
}


// Generates XML section related to the definition of a new goal
void
generateXMLSetNewGoal (FILE * outputFileXML, int node, double time,
			   double newGoalX, double newGoalY, double velocity,
			   double sideLength_WholeAreaX,
			   double sideLength_WholeAreaY)
{
	fprintf (outputFileXML, "\t\t<position_change>\n");
	fprintf (outputFileXML, "\t\t\t<node_id>%d</node_id>\n", node);
	fprintf (outputFileXML, "\t\t\t<start_time>%f</start_time>\n", time);
	fprintf (outputFileXML, "\t\t\t<destination>\n");
	fprintf (outputFileXML, "\t\t\t\t<xpos>%f</xpos>\n",
		 newGoalX - sideLength_WholeAreaX / 2.0);
	fprintf (outputFileXML, "\t\t\t\t<ypos>%f</ypos>\n",
		 newGoalY - sideLength_WholeAreaY / 2.0);
	fprintf (outputFileXML, "\t\t\t</destination>\n");
	fprintf (outputFileXML, "\t\t\t<velocity>%f</velocity>\n", velocity);
	fprintf (outputFileXML, "\t\t</position_change>\n");
}


int
main (int argc, char *argv[])
{

//----------------------------
//Definitions of the constants
//----------------------------

//number of hosts
	int numHosts = 100;

//connection threshold
	double connectionThreshold = 0.1;

//length of the sides of the area
	double sideLength_WholeAreaX = 2000.0;
	double sideLength_WholeAreaY = 2000.0;

//radius of the transmission area
	double radius = 200.0;

//number of rows
	int numberOfRows = 20;

//number of columns
	int numberOfColumns = 20;

//simulationTime
	double totalSimulationTime = 1000.0;

//refreshTime
	double stepInterval = 1.0;

//low bound speed of the host
	double minHostSpeed = 1.0;
	double maxHostSpeed = 10.00;

//reconfiguration interval
//it defines the interval between two reconfigurations
	double reconfigurationInterval = 10000.0;

//rewiring probability
	double rewiringProb = 0.1;

//numberOfGroups
	int numberOfGroups = 10;

//simulation time
	double simTime = 0.0;

//number of travellers
	int numberOfTravellers = 0;

//speed of the travellers
	double travellerSpeed = 20.0;

// HIBOP: probability of remaining in the same group for a while
	double travProb = 0;

//verbose on/off
	bool verbose = false;

//generateXML on/off
	bool generateXML = false;

//girvan Newman clustering algorithm on off
	bool girvanNewmanOn = true;

//colocation traces on/off
	bool colocationTracesOn = true;

//communities traces on/off
	bool communitiesTracesOn = true;
	string adjacencyFile ("communities.idat");	// HIBOP

//deterministic on/off
	bool deterministicOn = true;
	bool nsTraceOn = false;

//drift
	float drift;

//debug variables
	bool firstTime = true;
	double timeFirstReconfiguration;

// input from file - HIBOP
	string inputFileName;
	bool inputFromFile = false;

	string workingFolder (".");

	// Probability of remaining in a non-home cell
	double remainingProb = 0;

//***********************
// Retrieving input values
// ***********************


	char ch;
	int inputSeed;
	vector < bool > getOpt (18, false);

	while ((ch = getopt (argc, argv, "hxNy:n:s:S:p:l:d:X:Y:c:C:R:G:r:T:t:w:v:g:a:d:A:b:i:e:f:o:")) != EOF) {
		switch (ch) {
			case 'h':
				//shows helps
				header ();
				usage ();
				exit (1);
			case 'x':
				//generate XML file using ANSim format
				generateXML = true;
				break;
			case 'n':
				//sets the number of hosts
				numHosts = atoi (optarg);
				getOpt[0] = true;
				break;
			case 's':
				//sets lower bound of the speed of the hosts
				minHostSpeed = (double) atof (optarg);
				getOpt[2] = true;
				break;
			case 'S':
				//sets upper bound of the speed of the hosts
				maxHostSpeed = (double) atof (optarg);
				getOpt[3] = true;
				break;
			case 'p':
				//sets the connection threshold
				connectionThreshold = (double) atof (optarg);
				getOpt[4] = true;
				break;
			case 'X':
				//sets the side length of the whole area (x coordinates)
				sideLength_WholeAreaX = (double) atof (optarg);
				getOpt[5] = true;
		
				break;
			case 'Y':
				//sets the side length of the whole area (y coordinates)
				sideLength_WholeAreaY = (double) atof (optarg);
				getOpt[6] = true;
				break;
			case 'R':
				//sets the number of rows
				numberOfRows = atoi (optarg);
				getOpt[7] = true;
				break;
			case 'C':
				//sets the number of columns
				numberOfColumns = atoi (optarg);
				getOpt[8] = true;
				break;
			case 'T':
				//sets the transmission range (radius)
				radius = (double) atof (optarg);
				getOpt[9] = true;
				break;
			case 't':
				//sets total simulation time
				totalSimulationTime = (double) atof (optarg);
				getOpt[1] = true;
				break;
			case 'r':
				//sets the reconfiguration interval
				reconfigurationInterval = (double) atof (optarg);
				getOpt[10] = true;
				break;
			case 'w':
				//sets the rewiring probability
				rewiringProb = (double) atof (optarg);
				getOpt[11] = true;
				break;
			case 'G':
				//sets the number of groups of the social network
				numberOfGroups = atoi (optarg);
				getOpt[13] = true;
				break;
			case 'g':
				//sets the seed of the random number generator
				inputSeed = atoi (optarg);
				getOpt[14] = true;
				setSeed (inputSeed);
				break;
			case 'c':
				//sets the number of carriers/travellers
				numberOfTravellers = atoi (optarg);
				getOpt[15] = true;
				break;
			case 'v':
				//sets the speed of the travellers
				travellerSpeed = (double) atof (optarg);
				getOpt[16] = true;
				break;
			case 'a':
				//sets the Girvan/Newman algorithm on/off
				if (strcmp (optarg, "on") == 0)
				girvanNewmanOn = true;
				else
				girvanNewmanOn = false;
				break;
			case 'A':
				//sets colocation traces on/off
				if (strcmp (optarg, "on") == 0)
				colocationTracesOn = true;
				else
				colocationTracesOn = false;
				break;
			case 'b':
				//sets communities on/off
				if (strcmp (optarg, "on") == 0)
				communitiesTracesOn = true;
				else
				communitiesTracesOn = false;
				break;
			case 'i':
				inputFileName = optarg;
				inputFromFile = true;
				break;
			case 'e':
				adjacencyFile = optarg;
				break;
			case 'f':
				workingFolder = optarg;
				break;
			case 'd':
				//sets deterministic on/off
				if (strcmp (optarg, "on") == 0)
					deterministicOn = true;
				else {
					cerr << ">> Probabilistic Selection of next goal <<\n";
					deterministicOn = false;
				}
				break;
			case 'o':
				// set probability of remaining in a non-home cell
				remainingProb = (double) atof (optarg);
				getOpt[17] = true;
				break;
			case 'N':
				nsTraceOn = true;
				break;
			case ':':
				//fprintf(stderr,"Option -%c requires an operand\n", optopt);
				usage ();
				exit (1);
			case '?':
				//
				usage ();
				exit (1);
			default:
				break;
		}
	}

	if (inputFromFile == true) {
		parsingConfiguration (numHosts, minHostSpeed, maxHostSpeed,
				  connectionThreshold, sideLength_WholeAreaX,
				  sideLength_WholeAreaY, numberOfRows,
				  numberOfColumns, radius, totalSimulationTime,
				  reconfigurationInterval, rewiringProb,
				  numberOfGroups, inputSeed, numberOfTravellers,
				  travellerSpeed, travProb, remainingProb, inputFileName,
				  getOpt);
	}


	// Set seed
	setSeed (inputSeed);

	bool isConnected[numHosts][numHosts];
	bool isConnectedWithAP[numHosts];

	hostsItem hosts[numHosts];
	hostsItem accessPoint; // ap element


	cellsItem cells[numberOfRows][numberOfColumns];

	double sideLengthX;
	double sideLengthY;

	//cell attractivity
	float CA[numberOfRows][numberOfColumns];

	//Variables used for the generation of the mobility scenario using the results
	//of complex networks theory and in particular the Girvan-Newman algorithm

	//interaction->adjacency threshold
	double threshold = connectionThreshold;

	int numberOfMembers[numHosts];
	int **adjacency;
	double **interaction;
	int **groups;

	//output stream
	std::ofstream out1;

	//input fileName
	char resultsFileName[] = "results.txt";

	//double last values registered
	double lastValues[numHosts][numHosts];
	double lastValuesWithAP[numHosts];

	//array [hostId,communitiesId]
	int communities[numHosts][2];

	//probability of moving to the cell [c][r]
	float a[numberOfRows][numberOfColumns];

	struct probRange {
		float min;
		float max;
	};

	probRange p[numberOfRows][numberOfColumns];

	//probability space distribution in the range [0..1]
	
	// --- CA
	int newComm; // The community the node is currently in
	vector< int > currComm(numHosts); // The community the node was in during last movement
	vector< double > currCommStartTime(numHosts); // Time the node has entered the community
	vector<vector <double> > frequencies(numHosts, vector<double>(numberOfRows*numberOfColumns, 0));
	vector<int> previousGoalCommunity(numHosts);
	vector<vector <bool> > friendCommunity(numHosts, vector<bool>(numberOfRows*numberOfColumns));
	// ---

	// **********************************   
	// Initialisation of the output files   
	// **********************************

#ifdef WIN32
	LPSECURITY_ATTRIBUTES attr;
	attr = NULL;
	if (workingFolder != ".") {
		CreateDirectory (workingFolder.c_str (), attr);
	}
	CreateDirectory ((workingFolder + "/idat").c_str (), attr);
#else
	// UNIX CODE
	mkdir (workingFolder.c_str (), S_IRWXU | S_IRWXG | S_IRWXO);
	mkdir ((workingFolder + "/idat").c_str (), S_IRWXU | S_IRWXG | S_IRWXO);
#endif

	//Creation of the ns-2 traces output file
	//if ((outputFile = fopen ("socMov.tr","wt"))==NULL)
	//		fprintf (stderr,"Cannot open %s\n","socMov.tr");
	
	// Output File for Omnet++
	fstream outputOmnet;
	if ( true) {
		stringstream fileName;
		fileName << workingFolder << "/contacts.txt";
		outputOmnet.open(fileName.str().c_str(), ios::out);
	}

	//Creation of the output file XML
	FILE *outputFileXML;
	generateXML = false;
	if (generateXML == true)
		if ((outputFileXML = fopen ("socMov.xml", "wt")) == NULL)
			fprintf (stderr, "Cannot open %s\n", "socMov.xml");


	//Creation of the colocation statistics files
	FILE *outputColocFiles[numHosts+1];
	if (colocationTracesOn == true)
		for (int i = 0; i < numHosts+1; i++) {
			stringstream fileName;
			fileName << workingFolder << "/idat/" << i << ".idat";
			if ((outputColocFiles[i] =  fopen (fileName.str ().c_str (), "wt")) == NULL)
				//fprintf(stderr,"Cannot open %s\n",fileName.str());
				cerr << "Cannot open " << fileName.str () << "\n";
		} //end for

	//Creation of the communities statistics files
	FILE *outputCommunitiesFile;
	if (communitiesTracesOn == true)
		if ((outputCommunitiesFile = fopen ((workingFolder + "/idat/" + adjacencyFile).c_str (), "wt")) == NULL)
			fprintf (stderr, "Cannot open communities.idat\n");

	// --- CA: Creation of the community awareness statistics files
	FILE *outputCommunityAwarenessFiles[numHosts];
	if (colocationTracesOn == true)
		for (int i = 0; i < numHosts; i++) {
			stringstream fileName;
			fileName << workingFolder << "/idat/" << i << "_comm.idat";
			if ((outputCommunityAwarenessFiles[i] =  fopen (fileName.str ().c_str (), "wt")) == NULL)
				//fprintf(stderr,"Cannot open %s\n",fileName.str());
				cerr << "Cannot open " << fileName.str () << "\n";
		} //end for

#if FREQUENTIAL_WEIGHTS
	FILE *outputFreqFiles[numHosts];
	for (int i = 0; i < numHosts; i++) {
		stringstream fileName;
		fileName << workingFolder << "/idat/" << i << "_freq.idat";
		if ((outputFreqFiles[i] =  fopen (fileName.str ().c_str (), "wt")) == NULL)
			//fprintf(stderr,"Cannot open %s\n",fileName.str());
			cerr << "Cannot open " << fileName.str () << "\n";
	} //end for
#endif

#if SOCIAL_WEIGHTS
	FILE *outputSocialFiles[numHosts];
	for (int i = 0; i < numHosts; i++) {
		stringstream fileName;
		fileName << workingFolder << "/idat/" << i << "_soc.idat";
		if ((outputSocialFiles[i] =  fopen (fileName.str ().c_str (), "wt")) == NULL){
			//fprintf(stderr,"Cannot open %s\n",fileName.str());
			cerr << "Cannot open " << fileName.str () << "\n";
		} else {
		//	cerr << "Social files created.\n";
		}
	} //end for
#endif

#if DISTANCE_DISTRIBUTION
	// OUTPUT FILE for distance distribution	
	FILE *outputDistanceDistr[numHosts];
	for (int i = 0; i < numHosts; i++) {
		stringstream fileName;
		fileName << workingFolder << "/idat/" << i << "_distDistr.idat";
		if ((outputDistanceDistr[i] =  fopen (fileName.str ().c_str (), "wt")) == NULL) {
			//fprintf(stderr,"Cannot open %s\n",fileName.str());
			cerr << "Cannot open " << fileName.str () << "\n";
		}
	} //end for
#endif
	
	// ***************************
	// Initialisation of the model
	//****************************
	sideLengthX = sideLength_WholeAreaX / ((double) numberOfRows);
	sideLengthY = sideLength_WholeAreaY / ((double) numberOfColumns);

	for (int i = 0; i < numberOfRows; i++)
		for (int j = 0; j < numberOfColumns; j++) {
			cells[i][j].minX = ((double) i) * sideLengthX;
			cells[i][j].minY = ((double) j) * sideLengthY;
			cells[i][j].numberOfHosts = 0;
		}

	//setup of the links
	for (int i = 0; i < numHosts; i++) {
		for (int l = 0; l < numHosts; l++) {
			isConnected[i][l] = false;
			lastValues[i][l] = 0.0;
			//nodeModule1->setIsInReach(l,false);
		}
		isConnectedWithAP[i] = false;
		lastValuesWithAP[i] = 0.0;
	}

	//initialization of the travellers
	for (int i = 0; i < numberOfTravellers; i++) {
		hosts[i].isATraveller = true;
		//definition of the initial speeds of the travellers			
		hosts[i].speed = travellerSpeed;
	}

	for (int i = numberOfTravellers; i < numHosts; i++) {
		hosts[i].isATraveller = false;

		// Use the formula proposed by Camp et a. for the stazionary version of RWP
                // Ref. Stationary distribution for the random waypoint mobility model
                double unif = generateRandomDouble(0, 1);
                hosts[i].speed = pow(maxHostSpeed, unif)/pow(minHostSpeed, unif-1);
		// ORIG_CMM: // hosts[i].speed = generateRandomDouble (minHostSpeed, maxHostSpeed);
	 }
	 cerr << "Camp's rule for stationary speed distribution DONE\n";

	accessPoint.speed = 0;

	double numberOfReconfigurations = 0.0;
	double nextReconfigurationTime = 0.0;

	int initialNumberOfGroups = numberOfGroups;

	// For saving the home cell information 
	int home_X[numberOfGroups];
	int home_Y[numberOfGroups];
	//Update of the positions

	//adjacency=initialise_int_array(numHosts);
	//interaction=initialise_double_array(numHosts);
	//groups=initialise_int_array(numHosts);

	for (simTime = 0.0; simTime < totalSimulationTime; simTime = simTime + stepInterval) {
		cerr << "Simulation Time : " << simTime << "\r";
		//reconfiguration mechanism
		if (simTime == nextReconfigurationTime) {

			for (int i = 0; i < numberOfRows; i++)
				for (int j = 0; j < numberOfColumns; j++) {
					cells[i][j].numberOfHosts = 0;
				}
	
			//start correction
			adjacency = initialise_int_array (numHosts);
			//print_int_array(adjacency,numHosts);
			interaction = initialise_double_array (numHosts);
	
			//print_double_array(interaction,numHosts);
	
			groups = initialise_int_array (numHosts);
			//end correction
	
			numberOfReconfigurations = numberOfReconfigurations + 1.0;
			nextReconfigurationTime = reconfigurationInterval * numberOfReconfigurations;
	
			bool splitted = true;
			double previousModth = 0.0;
			double modth = 0.1;
	
			if (girvanNewmanOn == true) {
	
				initialise_weight_array_ingroups (interaction, numHosts, initialNumberOfGroups, rewiringProb, threshold, inputSeed);
				generate_adjacency (interaction, adjacency, threshold, numHosts);
	
				//clustering using the Girvan-Newman algorithm
				do {
					for (int i = 0; i < numHosts; i++)
						numberOfMembers[i] = 0;
					splitted = false;
					double betw[numHosts];
					for (int i = 0; i < numHosts; i++)
						betw[i] = 0;
					calculate_betweenness (betw, adjacency, numHosts);
		
					for (int i = 0; i < numHosts; i++)
						numberOfMembers[i] = 0;
					numberOfGroups = getGroups (adjacency, groups, numberOfMembers, numHosts);
		
					if (verbose == true)
						printGroups (numberOfGroups, groups, numberOfMembers, numHosts);
					previousModth = modth;
		
					modth = splitNetwork_Threshold (adjacency, betw, numHosts, modth);
					if (verbose == true)
						cout << "\nModth is equal to " << modth;
		
				}
				while ((previousModth < modth) && (modth > -1));
				//print_double_array(interaction, numHosts);
			} else {			//end if (girvanNewmanOn==true)
			
				//communities based on the initial number of caves in the Caveman model		  
				//i.e., w=0
				initialise_weight_array_ingroups (interaction, numHosts, initialNumberOfGroups, 0, threshold, inputSeed);
				generate_adjacency (interaction, adjacency, threshold, numHosts);
	
				for (int i = 0; i < numHosts; i++)
					numberOfMembers[i] = 0;
				numberOfGroups = getGroups (adjacency, groups, numberOfMembers, numHosts);
				if (verbose == true)
					printGroups (numberOfGroups, groups, numberOfMembers, numHosts);
	
#if CONTROLLED_REWIRING
				initialise_weight_array_ingroups2(interaction, numHosts, initialNumberOfGroups, rewiringProb, threshold, inputSeed);
#else
				initialise_weight_array_ingroups(interaction, numHosts, initialNumberOfGroups, rewiringProb, threshold, inputSeed);

#endif
				generate_adjacency (interaction, adjacency, threshold, numHosts);
		
				//print_double_array(interaction, numHosts);
	
			} //end else

			int pointer = 0;
			vector < vector < bool > >usedCells (numberOfRows, vector < bool >  (numberOfColumns, false));
			vector < vector < int > >neighbors (numberOfRows, vector < int >(numberOfColumns, 0));
			int allowedNeigh = ALLOWED_NEIGH;
	
			for (int i = 0; i < numberOfGroups; i++) {
				/******** OUR CHANGE ********/
				// avoid that 2 groups are assigned to the same cell
				// avoid that 2 groups sare more than x edges
				int cellIdX;
				int cellIdY;
				bool neigh;
				do {
					neigh = false;
					cellIdX = generateRandomInteger (0, numberOfRows); // [1,4]
					cellIdY = generateRandomInteger (0, numberOfColumns);
					/*if ((cellIdX - 1 > 0) && ((usedCells[cellIdX - 2][cellIdY - 1] == true) 
						&& (neighbors[cellIdX - 2][cellIdY - 1] >= allowedNeigh))) {
						neigh = true;
					}
					if ((cellIdY - 1 > 0) && ((usedCells[cellIdX - 1][cellIdY - 2] == true)
						&& (neighbors[cellIdX - 1][cellIdY - 2] >= allowedNeigh))) {
						neigh = true;
					}
					if ((cellIdX - 1 < numberOfRows - 1) && ((usedCells[cellIdX][cellIdY - 1] == true)
						&& (neighbors[cellIdX][cellIdY - 1] >= allowedNeigh))) {
						neigh = true;
					}
					if ((cellIdY - 1 < numberOfColumns - 1) && ((usedCells[cellIdX - 1][cellIdY] == true)
						&& (neighbors[cellIdX - 1][cellIdY] >= allowedNeigh))) {
						neigh = true;
					}*/

					if (neighbors[cellIdX - 1][cellIdY - 1] >= allowedNeigh) { 
						neigh = true;
					}
				} while ( //(cellIdX == (numberOfRows/2+1) && cellIdY == (numberOfColumns/2+1)) ||
						(usedCells[cellIdX - 1][cellIdY - 1] == true || neigh == true));
	
				usedCells[cellIdX - 1][cellIdY - 1] = true;
				if (cellIdX - 1 > 0) 
					neighbors[cellIdX - 2][cellIdY - 1]++;
				if (cellIdY - 1 > 0)
					neighbors[cellIdX - 1][cellIdY - 2]++;
				if (cellIdX - 1 < numberOfRows - 1)
					neighbors[cellIdX][cellIdY - 1]++;
				if (cellIdY - 1 < numberOfColumns - 1)
					neighbors[cellIdX - 1][cellIdY]++;

				if (cellIdX - 1 > 0 && cellIdY-1 > 0)
					neighbors[cellIdX - 2][cellIdY - 2]++;
				if (cellIdX - 1 > 0 && cellIdY-1 < numberOfColumns -1)
					neighbors[cellIdX - 2][cellIdY]++;
				if (cellIdX - 1 < numberOfRows -1 &&  cellIdY-1 > 0)
					neighbors[cellIdX][cellIdY-2]++;
				if (cellIdX - 1 < numberOfRows -1 &&  cellIdY-1 < numberOfColumns - 1)
					neighbors[cellIdX][cellIdY]++;

				/******** END CHANGE ********/
	
				home_X[i] = cellIdX;
				home_Y[i] = cellIdY;
				 
	
				cerr << "The HOME cell for group" << i << " is " << (cellIdX - 1)*numberOfColumns + cellIdY - 1 << "\n";
	
				for (int j = 0; j < numberOfMembers[i]; j++) {
					int hostId = groups[i][j];
					hosts[hostId - 1].cellIdX = cellIdX;
					hosts[hostId - 1].cellIdY = cellIdY;
		
					communities[pointer][0] = hostId;
					communities[pointer][1] = i + 1;
					pointer++;
		
					//increment the number of the hosts in that cell
					cells[cellIdX - 1][cellIdY - 1].numberOfHosts += 1;
					
					// --- CA: Initialise current community
					currComm[hostId - 1] = (cellIdX-1) * numberOfColumns + (cellIdY-1);
					previousGoalCommunity[hostId - 1] = (cellIdX-1) * numberOfColumns + (cellIdY-1);
					currCommStartTime[hostId - 1]=simTime;
					//cerr <<  hostId - 1<< " " << currComm[hostId - 1] << " " << currCommStartTime[hostId - 1] << "\n";
					// ---
				}
			} // end for

			// Get friend communities for each node
				// Compute CA : being t=0, all nodes are in their home cell.
			for (int thisNode = 0; thisNode < numHosts; thisNode++){
				for (int c = 0; c < numberOfRows; c++)
					for (int r = 0; r < numberOfColumns; r++)
						CA[c][r] = 0;
		
				for (int n = 0; n < numHosts; n++) {
					int group_of_i = thisNode % numberOfGroups;
					int group_of_n = n % numberOfGroups;

					if ((n < thisNode) || (n > thisNode)) {
						
						CA[hosts[n].cellIdX - 1][hosts[n].cellIdY - 1] += adjacency[thisNode][n];
						
						/*
						// IF n and i belong to different groups 
						if ( group_of_n != group_of_i ) {
							// take into account the real "weight" of the external relationship
							// each node is accounted for when computing the attractivity of ITS home cell!
							CA[home_X[group_of_n] - 1][home_Y[group_of_n] - 1] += adjacency[thisNode][n];
						}
		
						// IF n and i belong to the same group, count n in the home's attractivity 
						if (group_of_n == group_of_i ) {
							//CA[home_X[group_of_n] - 1][home_Y[group_of_n] - 1] += adjacency[i][n];
							// take into account the affiliation to the same home community
							CA[home_X[group_of_i] - 1][home_Y[group_of_i] - 1] += 1;
						}*/
					}
				}
	
				double totAttractivity = 0;
				for (int c = 0; c < numberOfRows; c++) {
					for (int r = 0; r < numberOfColumns; r++) {
		
						if (cells[c][r].numberOfHosts != 0) {
							//CA[c][r] = CA[c][r] /	(double) cells[c][r].numberOfHosts;  //CA[c][r]=CA[c][r];	
							CA[c][r] = CA[c][r];
							friendCommunity[thisNode][c*numberOfColumns+r] = CA[c][r] > 0 ? true : false;

						} else {
							CA[c][r] = 0;
							friendCommunity[thisNode][c*numberOfColumns+r] = false;
						}

						totAttractivity += CA[c][r];

						if ( friendCommunity[thisNode][c*numberOfColumns+r] == true ) {
							if ( c*numberOfColumns+r ==
									(home_X[thisNode % numberOfGroups]-1) * numberOfColumns
									+ home_Y[ thisNode % numberOfGroups ]-1){
								/*cerr << "Il nodo " << thisNode << " ha relazione con la comunita' "
									<< c*numberOfColumns+r << " - Attractivity: " << CA[c][r] << " * \n";*/
							} else {
								/*cerr << "Il nodo " << thisNode << " ha relazione con la comunita' "
									<< c*numberOfColumns+r << " - Attractivity: " << CA[c][r] << "\n";*/

							}
						}
						//printf ("Deterministic %f ",CA[c][r]);
					}
				}
#if SOCIAL_WEIGHTS
				// Normalisation & output of weights
				for (int c = 0; c < numberOfRows; c++) {
					for (int r = 0; r < numberOfColumns; r++) {
						if ( c*numberOfColumns+r ==
									(home_X[thisNode % numberOfGroups]-1) * numberOfColumns
									+ home_Y[ thisNode % numberOfGroups ]-1){
							fprintf(outputSocialFiles[thisNode], "%d %f *\n", c*numberOfColumns+r, CA[c][r]/totAttractivity);
						} else {
							fprintf(outputSocialFiles[thisNode], "%d %f\n", c*numberOfColumns+r, CA[c][r]/totAttractivity);
						}
					}
				}
				//fclose(outputSocialFiles[thisNode]);
#endif
			}
			// -- END : Get friend communities for each node


			cerr << "End Initialisation\n";
	
			if (firstTime == true && communitiesTracesOn == true) {
	
				int temp1 = 0;
				int temp2 = 0;
				for (int i = 0; i < numHosts; i++) {
					for (int j = 0; j < numHosts - 1; j++) {
		
						if (communities[j][0] > communities[j + 1][0]) {
							temp1 = communities[j + 1][0];
							temp2 = communities[j + 1][1];
			
							communities[j + 1][0] = communities[j][0];
							communities[j + 1][1] = communities[j][1];
			
							communities[j][0] = temp1;
							communities[j][1] = temp2;
						}
					}		//end for (int j=0
				}		// end for (int i=0
				for (int i = 0; i < numHosts; i++) {
					fprintf (outputCommunitiesFile, "%d%s", communities[i][0], " ");
					fprintf (outputCommunitiesFile, "%d%s", communities[i][1], "\n");
				}
				fclose (outputCommunitiesFile);
			}			//end communitiesTracesOn
	
	
			if (firstTime == true) {
				cerr << "Camp's rule for stationary location distribution DONE\n";
				/* Insertion of an AP in the middle of the scenario. */
				accessPoint.currentX = sideLength_WholeAreaX/2;
				accessPoint.currentY = sideLength_WholeAreaY/2;
				// 2007/05/29 - End of changes

				//definition of the initial position of the hosts
				for (int k = 0; k < numHosts; k++) {
					/*
                                        // ORIG_CMM
                                        hosts[k].currentX = cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minX +
                                                                 generateRandomDouble (0.0, sideLengthX);
                                        hosts[k].currentY = cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minY +
                                                                 generateRandomDouble (0.0, sideLengthY);
                                        */

					// 26/11/07 - In the RWP version of CMM the initial positions must follow Camp's rule
					double r = 0;
                                        do {
                                                // x1 and y1
                                        hosts[k].currentX = cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minX +
                                                                 generateRandomDouble (0.0, sideLengthX);
                                        hosts[k].currentY = cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minY +
                                                                 generateRandomDouble (0.0, sideLengthY);

                                                // x2 and y2
                                        hosts[k].goalCurrentX = cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minX +
                                                                        generateRandomDouble (0.0, sideLengthX);
                                        hosts[k].goalCurrentY =  cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minY +
                                                                        generateRandomDouble (0.0, sideLengthY);

                                        r = pow(pow(hosts[k].goalCurrentX - hosts[k].currentX, 2) +
                                                        pow(hosts[k].goalCurrentY - hosts[k].currentY, 2), 1/2) / sqrt(2);

                                        } while (generateRandomDouble(0,1) >= r);
                                        double unif2 = generateRandomDouble(0,1);
                                        hosts[k].currentX = unif2 * hosts[k].currentX + (1-unif2) * hosts[k].goalCurrentX;
                                        hosts[k].currentY = unif2 * hosts[k].currentY + (1-unif2) * hosts[k].goalCurrentY;
                                        
                                        hosts[k].previousGoalX = hosts[k].currentX;
                                        hosts[k].previousGoalY = hosts[k].currentY;
                                        
                                        //cerr << "Initial assignment for node " << k << " DONE\n";

					// END - Camp's rule

				//	cerr << "Posizione Iniziale NODO " << k << ": [" << hosts[k].currentX
				//<< ", " << hosts[k].currentY << "]\n";
				}
				
				cerr << " ====== DONE FOR ALL =======\n";

	//=========================================================================================================================
				timeFirstReconfiguration = simTime;
				firstTime = false;
		
				// ************************
				// Generation of the traces
				// ************************
		
				if (nsTraceOn == true) {
					printf ("%s", "#Mobility Patterns Generator for ns-2 simulator\n");
					printf ("%s", "#Based on the Community Based Mobility Model\n");
					printf ("%s", "#Copyright (C) Mirco Musolesi University College London\n");
					printf ("%s", "#m.musolesi@cs.ucl.ac.uk\n");
					printf ("%s", "#Version 0.3 December 2006\n");
		
					printf ("%s", "set god_ [God instance]\n");
		
					//Generating initial positions of the hosts
					for (int i = 0; i < numHosts; i++) {
						printf ("%s%d%s%f%s", "$node_(", i, ") set X_ ", hosts[i].currentX, "\n");
						printf ("%s%d%s%f%s", "$node_(", i, ") set Y_ ", hosts[i].currentY, "\n");
						printf ("%s%d%s%f%s", "$node_(", i, ") set Z_ ", 0.0, "\n");
					}
					
					printf ("%s%d%s%f%s", "$node_(", numHosts, ") set X_ ", accessPoint.currentX, "\n");
					printf ("%s%d%s%f%s", "$node_(", numHosts, ") set Y_ ", accessPoint.currentY, "\n");
					printf ("%s%d%s%f%s", "$node_(", numHosts, ") set Z_ ", 0.0, "\n");

					printf ("%s %f %s%d%s %f %f %f%s", "$ns_ at", simTime, "\"$node_(", numHosts, ") setdest",
							accessPoint.currentX, accessPoint.currentY, 0.0, "\"\n");


				}
				//Generating initial positions of the hosts - XML code  
				if (generateXML == true) {
		
					generateXMLHeader (outputFileXML, sideLength_WholeAreaX, sideLength_WholeAreaY, radius, numHosts);
					fprintf (outputFileXML, "\t<node_settings>\n");
		
					for (int i = 0; i < numHosts; i++) {
						generateXMLSetNodePosition (outputFileXML, i, hosts[i].currentX, hosts[i].currentY, sideLength_WholeAreaX, sideLength_WholeAreaY);
					}
		
					fprintf (outputFileXML, "\t</node_settings>\n");
					fprintf (outputFileXML, "\t<mobility>\n");
		
				}		//end generation XML code
			}			//end if (firstTime==true)
	
			//definition of the initial goals
			if (firstTime == false) { // when firstTime == true this process follow the Camp's rule (see above)
				for (int k = 0; k < numHosts; k++) {
					hosts[k].goalCurrentX = cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minX +
									generateRandomDouble (0.0, sideLengthX);
					hosts[k].goalCurrentY =	 cells[hosts[k].cellIdX - 1][hosts[k].cellIdY - 1].minY +
									generateRandomDouble (0.0, sideLengthY);
				}
			}
	
	
			//generation of the traces - setting of the goals
			for (int i = 0; i < numHosts; i++) {
	
				hosts[i].absSpeed = hosts[i].speed;
	
				if (nsTraceOn == true) {
					printf ("%s %f %s%d%s %f %f %f%s", "$ns_ at", simTime, "\"$node_(", i, ") setdest",
							hosts[i].goalCurrentX, hosts[i].goalCurrentY, (hosts[i].absSpeed) / stepInterval, "\"\n");
		
				}
	
				if (generateXML == true)
					generateXMLSetNewGoal (outputFileXML, i, simTime, hosts[i].goalCurrentX, hosts[i].goalCurrentY,
						 hosts[i].absSpeed / stepInterval, sideLength_WholeAreaX, sideLength_WholeAreaY);
			} //end for (int i=0;i<numHosts;i++)
		} //end reconfiguration mechanism
	
		for (int i = 0; i < numHosts; i++) {
			// --- CA: Detect node's community - START
			// Check node's position
			newComm = (int)floor(hosts[i].currentX/sideLengthX) * numberOfColumns + (int)floor(hosts[i].currentY/sideLengthY);
			
			if (newComm != currComm[i]) {
				// ... if the node has changed its cell
				//cerr << i << " " << newComm << " " << hosts[i].currentX << " " << sideLengthX <<  " " <<
				//	(int)floor(hosts[i].currentX/sideLengthX) <<  " " << numberOfColumns << " " <<
				//	hosts[i].currentY << " " << sideLengthY << " " <<  (int)floor(hosts[i].currentY/sideLengthY) << "\n";
					// output: community start_time end_time
				//if ( currComm[i] == previousGoalCommunity[i] ) { 
				if (colocationTracesOn == true) {
					if ( friendCommunity[i][currComm[i]] == true ) { 
						fprintf (outputCommunityAwarenessFiles[i], "%d %d %d\n", currComm[i], (long) currCommStartTime[i],(long) simTime);
					} else {
						//fprintf (outputCommunityAwarenessFiles[i], "%d %d %d\n", currComm[i], (long) currCommStartTime[i],(long) simTime);
					}
				}

				frequencies[i][currComm[i]] += simTime - currCommStartTime[i];
				currComm[i] = newComm;
				currCommStartTime[i]=simTime;
			}
			// ---
	
			// Detect node's community - END
			
			if ((hosts[i].currentX > hosts[i].goalCurrentX + hosts[i].speed)
			     || (hosts[i].currentX <  hosts[i].goalCurrentX - hosts[i].speed)
			     || (hosts[i].currentY > hosts[i].goalCurrentY + hosts[i].speed)
			     || (hosts[i].currentY < hosts[i].goalCurrentY - hosts[i].speed)) {

				//move towards the goal		 
				if (hosts[i].currentX < (hosts[i].goalCurrentX - hosts[i].speed))
					hosts[i].currentX = hosts[i].currentX + hosts[i].speed;
				if (hosts[i].currentX > (hosts[i].goalCurrentX + hosts[i].speed))
					hosts[i].currentX = (hosts[i].currentX) - hosts[i].speed;
				if (hosts[i].currentY < (hosts[i].goalCurrentY - hosts[i].speed))
					hosts[i].currentY = (hosts[i].currentY) + hosts[i].speed;
				if (hosts[i].currentY > (hosts[i].goalCurrentY + hosts[i].speed))
					hosts[i].currentY = (hosts[i].currentY) - hosts[i].speed;

				
			} else {
	
				int selectedGoalCellX = 0;
				int selectedGoalCellY = 0;
				int previousGoalCellX = hosts[i].cellIdX;
				int previousGoalCellY = hosts[i].cellIdY;
	
			if (deterministicOn == true) {
	
				//Algorithm of the selection of the new cell
				for (int c = 0; c < numberOfRows; c++)
					for (int r = 0; r < numberOfColumns; r++)
						CA[c][r] = 0;
		
				for (int n = 0; n < numHosts; n++)
					if ((n < i) || (n > i))
						CA[hosts[n].cellIdX - 1][hosts[n].cellIdY - 1] += interaction[i][n];
	
				for (int c = 0; c < numberOfRows; c++)
					for (int r = 0; r < numberOfColumns; r++) {
		
						if (cells[c][r].numberOfHosts != 0)
							CA[c][r] = CA[c][r] / 	(double) cells[c][r].numberOfHosts;  //CA[c][r]=CA[c][r];	
						else
							CA[c][r] = 0;
						//printf ("Deterministic %f ",CA[c][r]);
					}
	
				int selectedGoalCellX2 = 0;
				int selectedGoalCellY2 = 0;
	
				double CAMax1 = 0;
				double CAMax2 = 0;
	
				for (int c = 0; c < numberOfRows; c++)
					for (int r = 0; r < numberOfColumns; r++)
						if (CA[c][r] > CAMax1) {
							//set the second best
							selectedGoalCellX2 = selectedGoalCellX;
							selectedGoalCellY2 = selectedGoalCellY;
							CAMax2 = CAMax1;
			
							selectedGoalCellX = c + 1;
							selectedGoalCellY = r + 1;
			
							CAMax1 = CA[c][r];
		
						} else if (CA[c][r] > CAMax2) {
							selectedGoalCellX2 = c + 1;
							selectedGoalCellY2 = r + 1;
							CAMax2 = CA[c][r];
						}
	
				/* ORIGINAL VERSION
				if ((previousGoalCellX==selectedGoalCellX)&&(previousGoalCellY==selectedGoalCellY)) {
					if (hosts[i].isATraveller==true) {
						if (selectedGoalCellX!=0) {
							selectedGoalCellX=selectedGoalCellX2;
							selectedGoalCellY=selectedGoalCellY2;
						}
					}
				}	*/
	
				// HIBOP - the traveller stays within a cell with probability p
				if (hosts[i].isATraveller == true) {
					if ((rand () / (RAND_MAX + 1.0)) < travProb) {
						selectedGoalCellX = previousGoalCellX;
						selectedGoalCellY = previousGoalCellY;
					} else {
						if ((previousGoalCellX == selectedGoalCellX) && (previousGoalCellY == selectedGoalCellY)) {
							if (selectedGoalCellX != 0) {
								selectedGoalCellX = selectedGoalCellX2;
								selectedGoalCellY = selectedGoalCellY2;
							}
						}
					}
				} // end if (hosts[i].isATraveller == true)
			}	/* end deterministic */ else {
				//probabilistic
				bool stampaDebug = false;
				if (i % 3 == 0) {
					stampaDebug = false;
					//stampaDebug = false;
				}
	
				/*****************************
				* Compute cell attractivity
				*****************************/
				if (!imInMyHome(previousGoalCellX, previousGoalCellY, home_X[i % numberOfGroups], home_Y[i % numberOfGroups])) {
					// If I am not in my home...
					if ((rand () / (RAND_MAX + 1.0)) < remainingProb) {
						// ... remain outside with probability *remainingProb*
						selectedGoalCellX = previousGoalCellX;
						selectedGoalCellY = previousGoalCellY;
					} else {
						// ... go back home with probability *remainingProb*
						selectedGoalCellX = home_X[i % numberOfGroups];
						selectedGoalCellY = home_Y[i % numberOfGroups];
					}
				} else { 
	
					//Algorithm of the selection of the new cell
					for (int c = 0; c < numberOfRows; c++)
						for (int r = 0; r < numberOfColumns; r++)
							CA[c][r] = 0.0;
	
					//if ( i == 0) cerr << "WRECOM07\n";
					// WRECOM 07 - BEGIN
					for (int n = 0; n < numHosts; n++) {
						int group_of_i = i % numberOfGroups;
						int group_of_n = n % numberOfGroups;

						if ((n < i) || (n > i)) {

							// IF n and i belong to different groups 
							if ( group_of_n != group_of_i ) {
								/*if ( i == 0 ) 
									cerr << i << " " << n << " DIFF GROUP - adj=" << adjacency[i][n] << "\n";*/

								// take into account the real "weight" of the external relationship
								// each node is accounted for when computing the attractivity of ITS home cell!
								CA[home_X[group_of_n] - 1][home_Y[group_of_n] - 1] += adjacency[i][n];
								
								//if ( i == 0 )
								//	cerr << "CA[" << (home_X[group_of_n] - 1)*numberOfColumns + home_Y[group_of_n] - 1
								//		<< "]=" << CA[hosts[n].cellIdX - 1][hosts[n].cellIdY - 1] << "\n";
							}
			
							// IF n and i belong to the same group, count n in the home's attractivity 
							if (group_of_n == group_of_i ) {
								/*if ( i == 0 ) 
									cerr << i << " " << n << " SAME GROUP - adj=" << adjacency[i][n] << "\n";*/

								CA[home_X[group_of_i] - 1][home_Y[group_of_i] - 1] += adjacency[i][n];
								
								// take into account the affiliation to the same home community
								//CA[home_X[group_of_i] - 1][home_Y[group_of_i] - 1] += 1;

								//if ( i == 0 )
								//	cerr << "CA[" << (home_X[group_of_i] - 1)*numberOfColumns + home_Y[group_of_i] - 1
								//		<< "]=" << CA[home_X[group_of_i] - 1][home_Y[group_of_i] - 1] << "\n";
							}
						}
						if (stampaDebug == true)
							cerr << adjacency[i][n] << " ";
					} // end for (int n = 0; n < numHosts; n++)
					// WRECOM 07 - END


					if (stampaDebug == true)
						cerr << "\nPrint CA for node " << i + 1 << ": \t\tSimTime=" << simTime << "\n";
	
					for (int c = 0; c < numberOfRows; c++) {
						for (int r = 0; r < numberOfColumns; r++) {
		
							if (cells[c][r].numberOfHosts != 0) {
								//CA[c][r]=CA[c][r]/(double)cells[c][r].numberOfHosts;
								CA[c][r] = CA[c][r];
								/*if ( i == 0 )
									cerr << "CA[" << c*numberOfColumns + r 	<< "]="	<< CA[c][r] << "\n";*/
							} else {
								CA[c][r] = 0;
							}
		
							if (stampaDebug == true) fprintf (stderr, "%f ", CA[c][r]);
		
						}
		
						if (stampaDebug == true) cerr << "\n";
					}
		
					/* Assegno ad ogni cella un intervallo appartente a [0,1] tanto più grande quanto più
					* grande la sua CA. L'unione di questi intervalli copre tutto [0,1]. Dopo
					* estraggo un numero random compreso tra 0 e 1 e seleziono come goal successivo
					* quella cella il cui intervallo associato contiene questo numero.
					*/
	
					//Denominator for the normalization of the values
					float denNorm = 0.00;
		
					// HIBOP
					//drift=0.01;   
					drift = 0.0;
		
		
					// if (stampaDebug == true) cerr << "denNorm \n";
		
					for (int c = 0; c < numberOfRows; c++) {
						for (int r = 0; r < numberOfColumns; r++) {
							denNorm = denNorm + CA[c][r] + drift;
							//if (stampaDebug == true)		cerr << denNorm << " ";
						}
						//if (stampaDebug == true) cerr << "\n";
					}
	
					if (stampaDebug == true) cerr << "denNorm = " << denNorm << "\n";
	
					for (int c = 0; c < numberOfRows; c++) {
						for (int r = 0; r < numberOfColumns; r++)
							if (CA[c][r] == 0)
								a[c][r] = drift / denNorm;
							else
								a[c][r] = (CA[c][r] + drift) / (+denNorm);
					}
	
					float current = 0.0;
					for (int c = 0; c < numberOfRows; c++){
						for (int r = 0; r < numberOfColumns; r++) {
							p[c][r].min = current;
			
							p[c][r].max = p[c][r].min + a[c][r];
							current = p[c][r].max;
		
						}
					}
	
					for (int c = 0; c < numberOfRows; c++){
						for (int r = 0; r < numberOfColumns; r++) {
		
							p[c][r].min = p[c][r].min;
							p[c][r].max = p[c][r].max;
			
							if (stampaDebug == true) 
								fprintf (stderr, "%s %d %d %s %f %f %s",
										 "Square", c, r, " : ", p[c][r].min, p[c][r].max, "\n");
						}
					}
	
					/* HIBOP - Debug */
					double CAMax1 = 0;
					int max_attract_cellX;
					int max_attract_cellY;
		
					for (int c = 0; c < numberOfRows; c++){
						for (int r = 0; r < numberOfColumns; r++)
							if (CA[c][r] > CAMax1) {
			
								max_attract_cellX = c + 1;
								max_attract_cellY = r + 1;
			
								CAMax1 = CA[c][r];
			
							}
					}
					/* Fine debug */
		
					//float infiniteDice= (float) random() / (float) 0x7fffffff;;
					float infiniteDice = (float) generateRandomDouble (0.0, 1.0);
		
					if (stampaDebug == true) fprintf (stderr, "%f%s", infiniteDice, "\n");
		
					for (int c = 0; c < numberOfRows; c++) {
						for (int r = 0; r < numberOfColumns; r++) {
							if ((infiniteDice > p[c][r].min) && (infiniteDice < p[c][r].max)) {
			
								selectedGoalCellX = c + 1;
								selectedGoalCellY = r + 1;
							}
						}
					}

					/* ===== Case 1 ======*/
		
					if (!imInMyHome(previousGoalCellX, previousGoalCellY, home_X[i % numberOfGroups], home_Y[i % numberOfGroups])) {
						// node wasn't in its home cell
						if (previousGoalCellX != selectedGoalCellX || previousGoalCellY != selectedGoalCellY) {
						// node has selected a goal in a different cell
							if ((rand () / (RAND_MAX + 1.0)) < remainingProb) {
								selectedGoalCellX = previousGoalCellX;
								selectedGoalCellY = previousGoalCellY;
							}
						}
					}

					/* ===== Case 2 =====*/
					// If the node was in its home-cell ...
					if (imInMyHome(previousGoalCellX, previousGoalCellY, home_X[i % numberOfGroups], home_Y[i % numberOfGroups])) {
						// ... and the new goal is outside ...
						if (!imInMyHome(selectedGoalCellX, selectedGoalCellY, 
									home_X[i % numberOfGroups], home_Y[i % numberOfGroups])) {
							// VVVV - debug - VVVV
							if (max_attract_cellX != home_X[i % numberOfGroups]
								&& max_attract_cellY != home_Y[i % numberOfGroups]) {
			
								// esco perché un'altra cella mi attrae di più
								//cerr << " ----------------------------  " << i + 1 <<
								//  " cambia cella perché un'altra lo attira di più !!!!! ------------------------\n";
							} else {
								// esco perché il meccanismo probabilistico ha selezionato un'altra cella
								//cerr << " ----------------------------  " << i + 1 <<
								//  " cambia cella a causa del meccanismo probabilistico !!!!! ------------------------\n";
							}
							// ^^^^ - debug - ^^^^
						}
					}

					/* ===== Case 3 ===== */
					// If the node was outside but it is now coming back to the home cell...
					if ((imInMyHome(selectedGoalCellX, selectedGoalCellY, home_X[i % numberOfGroups], home_Y[i % numberOfGroups]))
						&& (!imInMyHome(previousGoalCellX, previousGoalCellY, 
								home_X[i % numberOfGroups], home_Y[i % numberOfGroups]))) {
		
						// VVVVV Debug VVVVV
						//cerr << " ---------------------  " << i + 1 << " torna nella home !!!!! ---------------------\n";
		
					}
				}		// modifica 15/02	
			}		// end probabilistic
	
	
	
			//Re-definition of the number of hosts in each cell
			cells[previousGoalCellX - 1][previousGoalCellY - 1].numberOfHosts -= 1;
			cells[selectedGoalCellX - 1][selectedGoalCellY - 1].numberOfHosts += 1;

			previousGoalCommunity[i] = (previousGoalCellX-1) * numberOfColumns + (previousGoalCellY -1);
	
			double newGoalRelativeX = generateRandomDouble (0, sideLengthX);
			double newGoalRelativeY = generateRandomDouble (0, sideLengthY);
	
			//refresh of the information
			hosts[i].cellIdX = selectedGoalCellX;
			hosts[i].cellIdY = selectedGoalCellY;
		
#if DISTANCE_DISTRIBUTION
			int pc = (int)floor(hosts[i].previousGoalX/sideLengthX) * 
						numberOfColumns + (int)floor(hosts[i].previousGoalY/sideLengthY);
			int cc = (int)floor(hosts[i].goalCurrentX/sideLengthX) * 
						numberOfColumns + (int)floor(hosts[i].goalCurrentY/sideLengthY);
			double distance =  sqrt( pow(hosts[i].goalCurrentX - hosts[i].previousGoalX, 2) + 
									pow(hosts[i].goalCurrentY - hosts[i].previousGoalY, 2));
			
			
			/*cerr << "pc=" << pc << "[" << hosts[i].previousGoalX << "," << hosts[i].previousGoalY << "]\n";
			cerr << " cc=" << cc << "[" << hosts[i].goalCurrentX << "," << hosts[i].goalCurrentY << "]\n";
			cerr << " distance=" << distance << "\n";
			cerr << "=====================\n";*/
					

			fprintf(outputDistanceDistr[i], "%d %d %f\n", pc,cc, distance);
#endif
	
			hosts[i].previousGoalX = hosts[i].goalCurrentX;
			hosts[i].previousGoalY = hosts[i].goalCurrentY;
			
			hosts[i].goalCurrentX = cells[selectedGoalCellX - 1][selectedGoalCellY - 1].minX + newGoalRelativeX;
			hosts[i].goalCurrentY = cells[selectedGoalCellX - 1][selectedGoalCellY - 1].minY + newGoalRelativeY;
			hosts[i].absSpeed = hosts[i].speed;
	
			if (nsTraceOn == true) {
				printf ("%s %f %s%d%s %f %f %f%s", "$ns_ at", simTime, "\"$node_(", i, ") setdest", 
						hosts[i].goalCurrentX, hosts[i].goalCurrentY, (hosts[i].absSpeed) / stepInterval, "\"\n");
			}
	
			if (generateXML == true)
				generateXMLSetNewGoal (outputFileXML, i, simTime, hosts[i].goalCurrentX, hosts[i].goalCurrentY,
						hosts[i].absSpeed / stepInterval, sideLength_WholeAreaX, sideLength_WholeAreaY);
	
			}	
		}			// end for (int i=0;i<numHosts;i++) 
	
		for (int i = 0; i < numHosts; i++) {
			//update connectivity
			for (int j = 0; j < numHosts; j++) {
				//I am connected with me so the follwing condition must hold...:)
				if (i != j) {
		
					//calculation of the current distance
					double currentDistance = sqrt ((hosts[i].currentX - hosts[j].currentX) *
									(hosts[i].currentX - hosts[j].currentX) +
									(hosts[i].currentY - hosts[j].currentY) *
									(hosts[i].currentY - hosts[j].currentY));
		
					//if currentDistance<=radius then the hosts are connected
					if (currentDistance < radius) {
						//if the hosts has been previously disconnected, then they must be connected
						if (isConnected[i][j] == false) {
							isConnected[i][j] = true;
							lastValues[i][j] = simTime;
							
							outputOmnet << i << ":" << j << ":" << simTime << ":1\n";
		
						}
					} else {
						if (isConnected[i][j] == true) {
							if (simTime != 0) {
								//if the hosts has been previously connected, then they must be disconnected
								if (colocationTracesOn) {
									fprintf (outputColocFiles[i], "%d %d %d\n", j + 1,
											(long) lastValues[i][j],(long) simTime);
								}
								isConnected[i][j] = false;
								outputOmnet << i << ":" << j << ":" << simTime << ":0\n";
							}
						}
					} //end else 
				} //end if (i!=j)
			}
			/* 2007/05/29 - AP changes */
			// Do the same things for the AP!
			double currentDistance = sqrt((hosts[i].currentX - accessPoint.currentX) *
							(hosts[i].currentX - accessPoint.currentX) +
							(hosts[i].currentY - accessPoint.currentY) *
							(hosts[i].currentY - accessPoint.currentY));
			if (currentDistance < radius) {
				//if the host and the AP have been previously disconnected, then they must be connected
				if (isConnectedWithAP[i] == false) {
					isConnectedWithAP[i] = true;
					lastValuesWithAP[i] = simTime;
				}
			}  else { //else they are disconnected
				if (isConnectedWithAP[i] == true) {
					if (simTime != 0) {
					//if the hosts has been previously connected, then they must be disconnected
						if (colocationTracesOn) {
							fprintf(outputColocFiles[i], "%d %d %d\n", numHosts+1, (long) lastValuesWithAP[i],
									(long) simTime);
							fprintf(outputColocFiles[numHosts], "%d %d %d\n", i+1, (long) lastValuesWithAP[i],
									(long) simTime);
						}
						isConnectedWithAP[i] = false;
					}
				} 	// end if isConnectedWithAP[i][j] == true
			}	//end else
			
			/* 2007/05/29 - End AP changes */
		} // end connectivity for
	} // end simulation for

	/**********************
	 *  SIMULATION CLOSURE 
	 *  *******************/
	for (int i=0; i < numHosts; i++){

		// --- CA
		//if ( currComm[i] == previousGoalCommunity[i] || currComm[i] ==  (hosts[i].cellIdX-1)*numberOfColumns+(hosts[i].cellIdY-1)) { 
		if (colocationTracesOn == true) {
			if ( friendCommunity[i][currComm[i]] == true ) { 
				fprintf (outputCommunityAwarenessFiles[i], "%d %d %d\n", currComm[i], (long) currCommStartTime[i],(long) simTime);
			} else {
				//fprintf (outputCommunityAwarenessFiles[i], "%d %d %d\n", currComm[i], (long) currCommStartTime[i],(long) simTime);
			}
		}
		//fprintf (outputCommunityAwarenessFiles[i], "%d %d %d\n", currComm[i], (long) currCommStartTime[i],(long) simTime);
		//cerr << i << " " <<  currComm[i] << " " <<  (long) currCommStartTime[i] << " " << (long) simTime << "\n";
		//if (i == 0)
		//	cerr << previousGoalCommunity[i] << " " << (hosts[i].cellIdX-1)*numberOfColumns+(hosts[i].cellIdY-1) << "\n";
		frequencies[i][currComm[i]] += simTime - currCommStartTime[i];
		// ---

		for (int j = 0; j < numHosts; j++) {
			if (isConnected[i][j] == true) {
				if (simTime != 0) {
					//if the hosts has been previously connected, then they must be disconnected
					if (colocationTracesOn) {
						fprintf (outputColocFiles[i], "%d %d %d\n", j+1,(long) lastValues[i][j],(long) simTime);
					}
					isConnected[i][j] = false;
					outputOmnet << i << ":" << j << ":" << simTime << ":0\n";
				}
			}
		}

		if (isConnectedWithAP[i] == true) {
			//if the hosts has been previously connected, then they must be disconnected
			if (colocationTracesOn) {
				fprintf(outputColocFiles[i], "%d %d %d\n", numHosts+1, (long) lastValuesWithAP[i], (long) simTime);
				fprintf(outputColocFiles[numHosts], "%d %d %d\n", i+1, (long) lastValuesWithAP[i], (long) simTime);
			}
			isConnectedWithAP[i] = false;
		}
	}

#if FREQUENTIAL_WEIGHTS	
	// With this approach, weights are assigned based on the time spent in each community with which the node has a relationship
	// Check this part, I'm not sure results are consistent.
	for (int i=0; i < numHosts; i++){
		for (int j=0; j < numberOfRows * numberOfColumns; j++){
			// If the node has a social relation with this community, consider the time spent there
			if ( friendCommunity[i][j] == true ) { 
				if ( j ==  (home_X[i % numberOfGroups]-1) * numberOfColumns + home_Y[i % numberOfGroups]-1){
					fprintf(outputFreqFiles[i], "%d %f *\n", j, frequencies[i][j]/totalSimulationTime);
				} else {
					fprintf(outputFreqFiles[i], "%d %f\n", j, frequencies[i][j]/totalSimulationTime);
				}
			} else {
				// If the node has not any social relation, do not consider the time spend in the community
				fprintf(outputFreqFiles[i], "%d 0\n", j);
			}

		}
	}
#endif

	/* == FILE CLOSURE == */
	outputOmnet.close();
			
	for (int i = 0; i < numHosts; i++){
		
#if SOCIAL_WEIGHTS
		fclose(outputSocialFiles[i]);
#endif

#if FREQUENTIAL_WEIGHTS
		fclose(outputFreqFiles[i]);
#endif
		if ( colocationTracesOn ) {
			fclose(outputColocFiles[i]);
			fclose(outputCommunityAwarenessFiles[i]);
		}
		
#if DISTANCE_DISTRIBUTION
		fclose(outputDistanceDistr[i]);
#endif
	}
	//generation of the XML file - generation of the last lines
	if (generateXML == true) {
		fprintf (outputFileXML, "\t</mobility>\n");
		fprintf (outputFileXML, "\t<statistics>\n");
		fprintf (outputFileXML, "\t\t<stoptime>%f</stoptime>\n", simTime);
		fprintf (outputFileXML, "\t</statistics>\n");
		fprintf (outputFileXML, "</simulation>\n");
	}

	// close outputFileXML
	if (generateXML == true)
		fclose (outputFileXML);
	/* == END OF FILE CLOSURE == */
}

void
parsingConfiguration (int &numHosts, double &minHostSpeed,
			  double &maxHostSpeed,
			  double &connectionThreshold,
			  double &sideLength_WholeAreaX,
			  double &sideLength_WholeAreaY, int &numberOfRows,
			  int &numberOfColumns, double &radius,
			  double &totalSimulationTime,
			  double &reconfigurationInterval,
			  double &rewiringProb, int &numberOfGroups,
			  int &inputSeed, int &numberOfTravellers,
			  double &travellerSpeed, double &travProb, double &remainingProb,
			  string fileName, vector < bool > getOpt)
{
	vector < bool > missing (18, true);
	XMLNode xMainNode;
	XMLNode xNode;
	int n;
	int iteratorN = 0;
	int iteratorV = 0;
	string name;

// this open and parse the XML file:
	xMainNode = XMLNode::openFileHelper (fileName.c_str (), "HiBOp");

	xNode = xMainNode.getChildNode ("mobilityModel");

	n = xNode.nChildNode ("var");

	for (int i = 0; i < n; i++) {
	stringstream value;
	name = xNode.getChildNode ("var", &iteratorN).getAttribute ("name");
	value << xNode.getChildNode ("var",
					 &iteratorV).getAttribute ("value");
	if (name == "numNodes" && getOpt[0] == false) {
		value >> numHosts;
		missing[0] = false;
	} else if (name == "simulationLifeTime" && getOpt[1] == false) {
		value >> totalSimulationTime;
		missing[1] = false;
	} else if (name == "minHostSpeed" && getOpt[2] == false) {
		value >> minHostSpeed;
		missing[2] = false;
	} else if (name == "maxHostSpeed" && getOpt[3] == false) {
		value >> maxHostSpeed;
		missing[3] = false;
	} else if (name == "connectionThreshold" && getOpt[4] == false) {
		value >> connectionThreshold;
		missing[4] = false;
	} else if (name == "sideLength_WholeAreaX" && getOpt[5] == false) {
		value >> sideLength_WholeAreaX;
		missing[5] = false;
	} else if (name == "sideLength_WholeAreaY" && getOpt[6] == false) {
		value >> sideLength_WholeAreaY;
		missing[6] = false;
	} else if (name == "numberOfRows" && getOpt[7] == false) {
		value >> numberOfRows;
		missing[7] = false;
	} else if (name == "numberOfColumns" && getOpt[8] == false) {
		value >> numberOfColumns;
		missing[8] = false;
	} else if (name == "radius" && getOpt[9] == false) {
		value >> radius;
		missing[9] = false;
	} else if (name == "reconfigurationInterval" && getOpt[10] == false) {
		value >> reconfigurationInterval;
		missing[10] = false;
	} else if (name == "rewiringProb" && getOpt[11] == false) {
		value >> rewiringProb;
		missing[11] = false;
	} else if (name == "remainingProb" && getOpt[12] == false) {
		value >> remainingProb;
		missing[12] = false;
	} else if (name == "numberOfGroups" && getOpt[13] == false) {
		value >> numberOfGroups;
		missing[13] = false;
	} else if (name == "inputSeed" && getOpt[14] == false && inputSeed == 0) {
		cerr << "seed\n";
		value >> inputSeed;
		missing[14] = false;
	} else if (name == "numberOfTravellers" && getOpt[15] == false) {
		value >> numberOfTravellers;
		missing[15] = false;
	} else if (name == "travellerSpeed" && getOpt[16] == false) {
		value >> travellerSpeed;
		missing[16] = false;
	} else if (name == "travellerProb" && getOpt[17] == false) {
		value >> travProb;
		missing[17] = false;
	}
	}

	for (int i = 0; i < (int) missing.size (); i++) {
	if (missing[i] == true && getOpt[i] == false) {
		cerr << "CMM : Missing Arguments!!!" << i << "\n";
		exit (1);
	}
	}
}


bool
imInMyHome (int nowX, int nowY, int homeX, int homeY)
{
	return (nowX == homeX && nowY == homeY) ? true : false;
}
