/*

Library for analysis of social networks
	
Copyright (C) Mirco Musolesi University College London
              m.musolesi@cs.ucl.ac.uk

Version 0.1

May 2006

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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

The contact details of the author are the following:

Mirco Musolesi
Department of Computer Science - University College London
Gower Street London WC1E 6BT United Kingdom
Email: m.musolesi@cs.ucl.ac.uk
Phone: +44 20 7679 0391 Fax: +44 20 7387 1397
Web: http://www.cs.ucl.ac.uk/staff/m.musolesi


###
*/
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "socialNet.h"
using namespace std;

#define SEVEN_TRAV 0

const int array_sz = 1000;
const int infinity = -1;

ofstream fout;

// Simple print routines
//
void print_int_array(int **array, int array_size)
{	for (int i = 0; i < array_size; i++)
	{ cerr << array[i][0];
	  for (int j = 1; j < array_size; j++)
		cerr << ", " << array[i][j];
	  cerr << "\n";
	}
	cerr << "\n"; cerr.flush();
}

void print_double_array(double **array, int array_size)
{	for (int i = 0; i < array_size; i++)
	{ cerr << array[i][0];
	  for (int j = 1; j < array_size; j++)
		cerr << ", " << array[i][j];
	  cerr << "\n";
	}
	cerr << "\n"; cerr.flush();
}

void print_float_array(float **array, int array_size)
{	for (int i = 0; i < array_size; i++)
	{ cerr << array[i][0];
	  for (int j = 1; j < array_size; j++)
		cerr << ", " << array[i][j];
	  cerr << "\n";
	}
	cerr << "\n"; cerr.flush();
}


int **initialise_int_array(int array_size) {
	
	int **result = new int * [array_size];
	
	for (int i = 0; i < array_size; i++)
		result[i] = new int [array_size];

	return result;
}


double **initialise_double_array(int array_size)
{	double **result = new double * [array_size];
	
	for (int i = 0; i < array_size; i++)
		result[i] = new double [array_size];

	return result;
}



// This needs to be symmetric if there's no directionality to interactions
// Also, every node is considered to be adjacent to itself
// For the WS model, we start with a regular lattice, with n nodes, each of which is connected
// to its nearest k neighbours (k better be even). Then we rewire with probability p, given
// constraints that no two nodes are have more than one link between then, and a node is not linked
// to itself.
//
void initialise_adjacency_array_WS(int **adjacency, int n, int k, double p)
{	int i, j;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			if (i == j)
				adjacency[i][j] = 1;
			else	
				adjacency[i][j] = 0;

	// Create the regular lattice -- note that this is a circular lattice
	// If we iterate over the whole of the arrayand put in the k/2 forward links, then,
	// because this is symmetric, we will end up with the nearest k links in.
	//
	for (i = 0; i < n; i++)
		for (j = 1; j <= k/2; j++)
		{   int idx = (i+j) % n;		// Make it circular

			adjacency[i][idx] = 1;
			adjacency[idx][i] = 1;
		}

}



// This needs to be symmetric if there's no directionality to interactions
// Also, every node is considered to be adjacent to itself
//
void initialise_adjacency_array_random(int **adjacency, int array_size)
{	for (int i = 0; i < array_size; i++)
		for (int j = 0; j < array_size; j++)
			if (i == j)
				adjacency[i][j] = 1;
			else if (i < j)
			{ adjacency[i][j] = (((double) rand())/RAND_MAX < 0.5);
			  adjacency[j][i] = adjacency[i][j];
			}
}

void initialise_adjacency_array(int **adjacency, int array_size, int k, double p)
{	
	// initialise_adjacency_array_random(adjacency, array_size);
	initialise_adjacency_array_WS(adjacency, array_size, k, p);

}

// This needs to be symmetric if there's no directionality to interactions
// Also, every node is considered to be adjacent to itself so has weight zero
// Make this simple in the beginning, as model 1 -- l_i,j = 1 where i != j
void initialise_weight_array_WS(double **weight, int array_size)
{
	for (int i = 0; i < array_size; i++)
		for (int j = 0; j < array_size; j++)
			if (i == j)
				weight[i][j] = 0;
			else
				weight[i][j] = 1;
}

// This needs to be symmetric if there's no directionality to interactions
// Also, every node is considered to be adjacent to itself so has weight zero
//
void initialise_weight_array_random(double **weight, int array_size)
{	for (int i = 0; i < array_size; i++)
		for (int j = 0; j < array_size; j++)
			if (i == j)
				weight[i][j] = 0;
			else if (i < j)
			{ weight[i][j] = ((double) rand())/RAND_MAX;

			  // Impose arbitrary lower limit to ensure that reciprocals don't get
			  // unmanageably large
			  //
			  while (weight[i][j] < 0.01)
				weight[i][j] = ((double) rand())/RAND_MAX;

			  weight[j][i] = weight[i][j];
			}
}

//This class initialise the matrix creating a matrix composed of n disjointed groups with n rewiring
//threshold is the value under which the relationship is not considered important
void initialise_weight_array_ingroups(double ** weight, int array_size, int numberOfGroups, double probRewiring, double threshold, double seed) {

	//cerr << "Probabilità Rewiring = " << probRewiring << "\n";	
	srand((int)seed);
	
	int **groups=initialise_int_array(array_size);
	int numberOfMembers[numberOfGroups];
	
	for (int i=0;i<numberOfGroups;i++) {
		numberOfMembers[i]=0;
	}
	
	for (int i=0;i<array_size;i++) {
		int groupId=i%numberOfGroups;
		//cout<<"\n";
		//cout<<"Group id of host "<<i<<" is "<<groupId<<".\n";
		groups[groupId][numberOfMembers[groupId]]=i+1;
		numberOfMembers[groupId]+=1;	
	}

	int countRew=0;
        int countLinks=0;
        //cerr << "Prob Rew = " << probRewiring << "\n";
	
	int countPart=0;
	int countOut=0;                	

	bool done=false;	
	/* N.B. Weight matrix after rewiring is NOT symmetric */
	for (int i=0;i<array_size;i++){
		countPart=0;
		/*if (i%3 != 0){
			probRewiring=0;
		}else {
			//cerr << "Per il nodo " << i << "non azzero\n";
		}*/

		//if (done == false){ // used to force one node rewiring
		for (int j=0;j<array_size;j++) {
			if (i!=j ){
				if (areInTheSameGroup (i+1,j+1,groups,numberOfGroups,numberOfMembers)==true) {
					 countLinks++;	
					 if ((rand()/(RAND_MAX+1.0))<probRewiring) {	//rewiring
						done=true;
						countRew++;
						countPart++;
						bool found=false;
						for (int z=0;z<array_size;z++){
							if ((areInTheSameGroup (i+1,z+1,groups,numberOfGroups,numberOfMembers)==false)
								&&(weight[i][z]<threshold)&&(found==false)) {
								
								weight[i][z]=1.0-(rand()*threshold/(RAND_MAX+1.0));
								found=true;
								
							}
						}
						weight[i][j]=(rand()*threshold/(RAND_MAX+1.0));	
						//weight[i][j]=0;	
					} else {//no rewiring 
						weight[i][j]=1.0-(rand()*threshold/(RAND_MAX+1.0));
					}
			 	} else { //the hosts are not in the same cluster
			 
			 		weight[i][j]=(rand()*threshold/(RAND_MAX+1.0));
			 		//weight[i][j]=0;
			 	}
			}
		}
		//}
		/*if (countPart > 0){
			cerr << "Nodo " << i+1 << " Rewired " << countPart << "\n";
		}
		
		if (countPart >=array_size/numberOfGroups-2){
                        countOut++;
                }*/
	}

	//print_double_array(weight, array_size);

	/*if (countRew != 0) {
	                //cerr << "Tot rewired: " << countRew << "\t" << "Tot link: " << countLinks;
			cerr << countRew << " " << countLinks << "\n";
			cerr << ">= n-2 -> " << countOut << "\n";
	}*/                                        	
		
	
	//deletion of the variables
 	for (int i=0;i<array_size;i++) 
 		delete groups[i];
 	delete groups;

}

//This class initialise the matrix creating a matrix composed of n disjointed groups with n rewiring
//threshold is the value under which the relationship is not considered important
void initialise_weight_array_ingroups2(double ** weight, int array_size, int numberOfGroups, double probRewiring, double threshold, double seed) {

	cerr << "Probabilità Rewiring = " << probRewiring << "\n";	
	srand((int)seed);
	
	int **groups=initialise_int_array(array_size);
	int numberOfMembers[numberOfGroups];
	
	for (int i=0;i<numberOfGroups;i++) {
		numberOfMembers[i]=0;
	}
	
	for (int i=0;i<array_size;i++) {
		int groupId=i%numberOfGroups;
		//cout<<"\n";
		//cout<<"Group id of host "<<i<<" is "<<groupId<<".\n";
		groups[groupId][numberOfMembers[groupId]]=i+1;
		numberOfMembers[groupId]+=1;	
	}

	int countRew=0;
        int countLinks=0;
        //cerr << "Prob Rew = " << probRewiring << "\n";
	
	int countPart=0;
	int countOut=0;                	

	bool done=false;
		
	if (array_size/numberOfGroups < numberOfGroups ) {
		cerr << "ERRRRRRR\n";
		exit(99);
	}

	int trav=0;
	for (int gr = 1; gr < numberOfGroups; gr++){
		// Trav e gr are not from the same group
#if SEVEN_TRAV
		int num_trav_per_group = ( gr == 1) ? 4 : 3;
		for (int num_iter = 0; num_iter < num_trav_per_group; num_iter++){
#endif
			bool found=false;
			for (int z=0;z<array_size;z++){
				// find a node from the same group as trav ...
				if ((areInTheSameGroup (trav+1,z+1,groups,numberOfGroups,numberOfMembers)==true)
					&& (weight[trav][z] > threshold) && (trav != z) && (found==false) ) {
					
					// ... and reset their association
					weight[trav][z]=(rand()*threshold/(RAND_MAX+1.0));
					cerr << "Dissocio " << trav << " e " << z << " ... ";
					found=true;
				}
			}
			// Associate trav and gr
			weight[trav][gr]=1-(rand()*threshold/(RAND_MAX+1.0));
			cerr << "associo " << trav << " e " << gr << "\n";
			trav=trav+numberOfGroups;
#if SEVEN_TRAV
		}
#endif
	}
		
	//print_double_array(weight, array_size);

	/*if (countRew != 0) {
	                //cerr << "Tot rewired: " << countRew << "\t" << "Tot link: " << countLinks;
			cerr << countRew << " " << countLinks << "\n";
			cerr << ">= n-2 -> " << countOut << "\n";
	}*/                                        	
		
	
	//deletion of the variables
 	for (int i=0;i<array_size;i++) 
 		delete groups[i];
 	delete groups;

}



void rewire_WS(int **adjacency, double **weight, int n, int k, double p)
{
	//
	// OK, now rewire.
	// Start with the edges closest to each node (in a clockwise direction) and rewire with
	// probability p to give a new endpoint, provided that (i) it is not wired to us and
	// (ii) there is not already a link between us and that node. Then move out to the edges
	// connecting us to the next furthest away nodes (again in a clockwise direction) and repeat
	//
	for (int i = 1; i <= k/2; i++)
		for (int j = 0; j < n ; j++)
			if (((double) rand())/RAND_MAX <= p)
			{	int idx = (i+j) % n;	// Make it circular

				// Choose new endpoint
				//
				int new_end = rand() % n;

				while (new_end == idx || adjacency[j][new_end] != 0)
					new_end = rand() % n;

				adjacency[j][idx]	  = 0;
				adjacency[idx][j]	  = 0;
				adjacency[j][new_end] = 1;
				adjacency[new_end][j] = 1;
			}
}

void rewire_WS2(int **adjacency, double **weight, int n, int k, double p)
{
	//
	// OK, now rewire.
	// Start with the edges closest to each node (in a clockwise direction) and rewire with
	// probability p to give a new endpoint, provided that (i) it is not wired to us and
	// (ii) there is not already a link between us and that node. Then move out to the edges
	// connecting us to the next furthest away nodes (again in a clockwise direction) and repeat
	//
	for (int i = 1; i <= k/2; i++)
		for (int j = 0; j < n ; j++)
			if (((double) rand())/RAND_MAX <= p)
			{	int idx = (i+j) % n;	// Make it circular

				// Choose new endpoint
				//
				int new_end = rand() % n;

				while (new_end == idx || adjacency[j][new_end] != 0)
					new_end = rand() % n;

				adjacency[j][idx]	  = 0;
				adjacency[idx][j]	  = 0;
				adjacency[j][new_end] = 1;
				adjacency[new_end][j] = 1;

				weight[j][new_end]	  = 3;
				weight[new_end][j]    = 3;
			}
}



void initialise_weight_array(double **weight, int array_size)
{ 
	// initialise_weight_array_random(weight, array_size);
	initialise_weight_array_WS(weight, array_size);
}


void initialise_weight_array2(double ** weight, int array_size) {
		int seed = 5;
    	srand(seed);
		for (int i=0;i<array_size;i++)
			for (int j=0;j<array_size;j++) {
				if (i==j) weight[i][j]=1;
				if (i<j){
					weight[i][j]=rand()/(RAND_MAX+1.0);
					weight[j][i]=weight[i][j];
				}
			}
}


void rewire(int **adjacency, double **weight, int array_size, int k, double p)
{
	rewire_WS2(adjacency, weight, array_size, k, p);
}


// All pairs shortest paths problem.
//

void calculate_shortest_paths(double **shortest_path, double **shortest_path_ideal, int **adjacency, double **weight, int array_size)
{	int i, j, k;

	// ??? MIGHT BE BETTER TO USE DIJKSTRA HERE SINCE IT IS PROBABLE THAT THE NUMBER OF
	// EDGES << N^2
	
	// Use Floyd's algorithm
	//
	for (i = 0; i < array_size; i++)
		for (j = 0; j < array_size; j++)
		{	shortest_path[i][j]       = (adjacency[i][j] != 0 ? weight[i][j] : infinity);
			shortest_path_ideal[i][j] = weight[i][j];
		}	
		
	for (k = 0; k < array_size; k++)
		for (i = 0; i < array_size;i++)
		{	if (i == k)
				continue;

			for (j = 0; j < array_size; j++)
			{   if (j == i || j == k)
				  continue;


			  if (shortest_path[i][k] != infinity && shortest_path[k][j] != infinity)
			  { if (shortest_path[i][j] == infinity || (shortest_path[i][k] + shortest_path[k][j] < shortest_path[i][j]))
					shortest_path[i][j] = shortest_path[i][k] + shortest_path[k][j];

				if (shortest_path_ideal[i][j] == infinity || (shortest_path_ideal[i][k] + shortest_path_ideal[k][j] < shortest_path_ideal[i][j]))
					shortest_path_ideal[i][j] = shortest_path_ideal[i][k] + shortest_path_ideal[k][j];
			  }

			}

		}
}

// This routine now returns array_size * (array_size - 1) times the energy
// This is because this factor will be divided out in calculating E_glob in 
// any case.
//
double calculate_energy(double **shortest_path, int array_size)
{
  double sum = 0.0;

  for (int i = 0; i < array_size; i++)
	  for (int j = 0; j < array_size; j++)
		  if (i != j && shortest_path[i][j] != infinity)
		    sum += 1.0/shortest_path[i][j];

  // Whilst, strictly, we need this to calculate the energy
  // because it is normalised, we can remove it.
  //
  // sum /= array_size * (array_size - 1);
  return(sum);
}

// This routine now returns array_size * (array_size - 1) times the energy
// This is because this factor will be divided out in calculating E_loc in 
// any case.
//
// This energy is calculated in terms of local clustering -- take the immediate neighbours of a node
// and sum the reciprocal of the least cost paths between that little cluster. Then normalise by
// considering the case where all nodes are connected.
//
double calculate_energy_local(int node, int **adjacency, double **weight, int array_size)
{
	double res = 0.0;
	int	   idx_array_size = 0;

	static int    *idx_array					= new int[array_size];
	static int    **temp_adjacency				= initialise_int_array(array_size);
	static double **temp_weight					= initialise_double_array(array_size);
	static double **temp_shortest_path			= initialise_double_array(array_size);
	static double **temp_shortest_path_ideal	= initialise_double_array(array_size);

	// Calculate local energy for the graph formed by this node's immediate neighbours. Note that
	// this graph specifically excludes the node itself so as to measure fault tolerance -- how
	// efficient communication amongst the node's neighbours is if the node fails.
	//

	// work out the number of neighbours and their positions
	// the idx array can then be used to indirect through for further calcaulations
	//
	for (int i = 0; i < array_size; i++)
		if (adjacency[node][i] != 0 && i != node)		// i is an immediate neighbour of node
			idx_array[idx_array_size++] = i;

	// Fill in temporary arrays
	// Only need to fill in the entries we're interested in, which are in idx array.
	//
	for (int i = 0; i < idx_array_size; i++)
		for (int j = 0; j < idx_array_size; j++)
		{	int idx1 = idx_array[i];
			int idx2 = idx_array[j];

			temp_weight[i][j]	 = weight[idx1][idx2];
			temp_adjacency[i][j] = adjacency[idx1][idx2];
		}

	// Using only these nodes calculate shortest paths
	//
	calculate_shortest_paths(temp_shortest_path, temp_shortest_path_ideal, temp_adjacency, temp_weight, idx_array_size);

	// And work out the normalised global energy for this restricted set
	//
	double energy		= calculate_energy(temp_shortest_path, idx_array_size) ;
	double energy_ideal = calculate_energy(temp_shortest_path_ideal, idx_array_size);

	if (energy_ideal != 0)
		return (energy/energy_ideal);
	else
		return (0.0);
}

// E_glob is defined as E(G)/E_ideal(G)
//
double calculate_E_glob(double **shortest_path, double **shortest_path_ideal, int array_size)
{
	return (calculate_energy(shortest_path, array_size)/calculate_energy(shortest_path_ideal, array_size));
}	

// E_loc is defined as 1/N * sum(for all i) { E(G_i)/E_ideal(G_i) }
//
double calculate_E_loc(int **adjacency, double **weight, int array_size)
{	double sum = 0.0;

    for (int i = 0; i < array_size; i++)
		sum += calculate_energy_local(i, adjacency, weight, array_size);

	return (sum/array_size);
}	


// This is the cost evaluator function. As in the paper, it's the identity function for now.
//
double gamma(double l)
{ return l; }

// Cost is a very similar function to energy but uses weights instead of least cost paths
//
double calculate_cost(int **adjacency, double **weight, int array_size)
{	double numerator = 0.0;
	double denominator = 0.0;

	for (int i = 0; i < array_size; i++)
	  for (int j = 0; j < array_size; j++)
	  {	  if (i != j)
			{ denominator += gamma(weight[i][j]);
			  if (adjacency[i][j] != 0)
		        numerator += gamma(weight[i][j]);
			}
	  }

	  return (numerator/denominator);
}


//load matrix with double values with size array_size from a file called
//filename
int load_double_array(double** array,int array_size, const char* filename) {
	
	//no checks
	int i,j=0;
	int r;
	double value;
	//int * valueint;
	//FILE *fopen (const char *filename, const char *mode)
	FILE* file1=fopen (filename,"r");
	
	//cout<<"Loading matrix from file "<<filename<<"...\n";
	
	for (int i=0; i<array_size; i++) {
		for (int j=0; j<array_size; j++) {
			r=fscanf (file1, "%lf", &value);
			array[i][j]=value;
		}
	}
	
	fclose (file1);
}


//load matrix with double values with a file called filename
//the array_size is read in the file
//return the array size of the matrix
int load_double_array(double** array, const char* filename) {
	
	//no checks
	int i,j=0;
	int r;
	double value;
	//int * valueint;
	//FILE *fopen (const char *filename, const char *mode)
	FILE* file1=fopen (filename,"r");
	
	int array_size;
	
	r=fscanf(file1,"%d",&array_size);
	
	//<<"Loading matrix from file "<<filename<<"...\n";
	
	for (int i=0; i<array_size; i++) {
		for (int j=0; j<array_size; j++) {
			r=fscanf (file1, "%lf", &value);
			array[i][j]=value;
		}
	}
	
	fclose (file1);
	return array_size;
}


//generate the adjacency matrix from the weight matrix of size array_size given a certain threshold
void generate_adjacency (double** weightMat, int** adjacencyMat, double threshold, int array_size) {
	
	for (int i=0; i<array_size; i++)  
		for (int j=0; j<array_size; j++) {
			if (weightMat[i][j]>threshold)
				adjacencyMat[i][j]=1; 
			else
				adjacencyMat[i][j]=0;
		}
}


//print square matrix of integer of size array_size
void print_array_int (int** array, int array_size) {
	
	for (int i=0; i<array_size; i++) {
		for (int j=0; j<array_size; j++) { 
			printf ("%i  ",array[i][j]);
		}
		cout<<"\n";
	}
}
	
	
//print square matrix of double of size array_size
void print_array_double (double** array, int array_size) {
	for (int i=0; i<array_size; i++) {
		for (int j=0; j<array_size; j++) { 
			printf ("%.2f  ",array[i][j]);
		}
		cout<<"\n";
	}
}

void assign_distance_not_working(int current, int previous, double* distance, int d, int target, int* assigned, int**adjacency, int **pred,int array_size) {
	 		if (current!=target) {
	 			for (int k=0;k<array_size; k++)
	 			if ((k!=current)&&(k!=previous)) {	 					 			
	 				if (adjacency[current][k]==1) {					
	 					if (assigned[k]==1) {
	 						distance[k]=d+1;
	 						assigned[k]=0;
	 					}	
	 					else if (distance[k]==d) {
	 						//cout<<"Assign distance "<<d<<"\n";
	 						//assign_distance(k,current,distance,d,target,assigned,adjacency,pred,array_size); 
	 					}
	  				}
	 			}//if ((k!=current))
 		}
}	

//assign distance (used by calculate_betweenness)
void assign_distance(double* distance, int d, int* assigned, int**adjacency, int **pred, int* predNum, int array_size) {
 			
 			for (int k=0;k<array_size; k++)
	 			if (distance[k]==d) {
					for (int r=0;r<array_size;r++)	  			
	 					if (adjacency[k][r]==1) 
	 						if (assigned[r]==1) {
	 							distance[r]=d+1;
	 							assigned[r]=0;
	 							bool found=false;
	 							for (int p=0;p<predNum[r]; p++) 
	 								if (pred[r][p]==k) found=true;
	 							if (found==false) {
	 								pred[r][predNum[r]]=k;
	 								predNum[r]=predNum[r]+1;
	 							}
	 						}	
	 			}//if ((k!=current))
 		
}	

//calculate betweenness for a single vertex (used by calculate_betweenness)
void calculate_b (double* betw, double* distance, int **pred, int* predNum, int array_size) {
	
	for (int dist=15;dist>-1;dist--) 
		
		for (int k=0; k<array_size; k++)
			if (distance[k]==dist) {
				for (int i=0; i<predNum[k]; i++) {
					betw[pred[k][i]]=betw[pred[k][i]]+betw[k]/predNum[k];
				}
		
			}
}


//calculate betweenness
void calculate_betweenness(double* result_betw, int** adjacency, int array_size) {
	
	int* b=new int[array_size];
	double* distance= new double[array_size];
	int* assigned= new int[array_size];
	int* predNum= new int[array_size];
	double* betw= new double[array_size];
	double* current_betw=new double[array_size];
	
	int d;
	for (int i=0; i<array_size;i++) { 
		assigned[i]=1;
		distance[i]=15;
		betw[i]=1;
		current_betw[i]=0;
		//cout<<distance[i]<<"\n";
	
	}
		
	double **dist= initialise_double_array(array_size);
	int **pred= initialise_int_array(array_size);
		
	for (int i=0; i<array_size; i++) {
			for (int s=0; s<array_size;s++) { 
						assigned[s]=1;
						distance[s]=16;
						predNum[s]=0;
						betw[s]=1;
			}
			//for (int j=0; j<array_size; j++) {
			assigned[i]=0;
			distance[i]=0;
			for (int d=0;d<15;d++) 
				assign_distance(distance,d,assigned,adjacency,pred,predNum,array_size);			 
			//for (int j=0; j<array_size;j++)
			//	cout<<"Distance of "<<i+1<<" from "<<j+1<<" is "<<distance[j]<<"\n";
				
			//for (int l=0; l<array_size; l++) {
					
				//cout<<"The predecessors of "<<l+1<<" are: "; 
					
				//for (int m=0; m<predNum[l]; m++) {
					//cout<<pred[l][m]+1<<" ";
				//}
				//cout<<"\n";
			//}
				
			//}//end for (int j...
				
		
			calculate_b (betw,distance, pred, predNum,array_size); 		
		
			for (int z=0;z<array_size;z++) {
				current_betw[z]=current_betw[z]+betw[z];
			}	
				
	 }//end for (int i<0; i<array_size; i++)

	for (int z=0;z<array_size;z++)		
	result_betw[z]=current_betw[z];		

	//delete data structures
	for (int i = 0; i < array_size; i++) {
		delete pred[i];
		delete dist[i];
	}
	delete pred;
	delete dist;
	delete [] b;

	delete [] distance;
	delete []assigned;
	delete []predNum;
	delete []betw;
	delete []current_betw;
}


//assignToAGroup is a function used by getGroups
int assignToAGroup (int current, int currentGroup, int** adjacency, int** groups, int* numberOfMembers,int array_size,bool* assigned) {
	
	for (int k=0;k<array_size;k++) 
		
		if (adjacency[current][k]==1) 
			if (current!=k)
				if (assigned[k]==false) {
			 		numberOfMembers[currentGroup-1]=numberOfMembers[currentGroup-1]+1;
					groups[currentGroup-1][numberOfMembers[currentGroup-1]-1]=k+1;
					assigned[k]=true;
					assignToAGroup(k,currentGroup,adjacency, groups, numberOfMembers, array_size, assigned);	
				}
}


//return the number of groups given an adjacency matrix
//in groups
int getGroups (int** adjacency, int** groups, int* numberOfMembers, int array_size) {
	
	int numberOfGroups=0;
	bool assigned[array_size];
	
	for (int i=0; i<array_size; i++) {
		assigned[i]=false;	
	} 

	for (int i=0; i<array_size; i++) {
		
		if (assigned[i]==false) {
			numberOfGroups++;
			numberOfMembers[numberOfGroups-1]=numberOfMembers[numberOfGroups-1]+1;
			groups[numberOfGroups-1][numberOfMembers[numberOfGroups-1]-1]=i+1;
			assigned[i]=true;
			assignToAGroup(i,numberOfGroups,adjacency, groups, numberOfMembers, array_size, assigned);	
		}			
	}
	return numberOfGroups;
}


//Given a group of nodes and a node, stored in the array group of size equal
//to numberOfMembers, returns true if the node is in the group or false otherwise
bool isInGroup(int node, int* group, int numberOfMembers) {
	
	bool result=false;
	
	for (int k=0;k<numberOfMembers;k++)	
		if (group[k]==node) {
			result=true;
			//break;
		}
 	return result;
}


//split the social network into communities removing some link considering the betweeness of 
//each node
//returns the Modularity of the "splitting"
double splitNetwork (int** adjacency, double *betw,int array_size) { 

	int best=0;
	//double bestBetw=0;
	
	//double betw2[array_size];
	
	//calculate highest betweenness
	//for (int i=0;i<array_size;i++) {
	//	if (betw[i]>bestBetw) { 
	//		bestBetw=betw[i];
	//		best=i;	
	//	}	
	//}
	
	//bool oneFound=false;
	//double minBetw=bestBetw;
	//int linkToDelete=0;
	//for (int i=0;i<array_size;i++) {
	//	if ((adjacency[best][i]==1)&&(best!=i)){
	//				adjacency[best][i]=0;
	//				adjacency[i][best]=0;	
	//				calculate_betweenness(betw2, adjacency, array_size);
	//				if (betw2[best]<minBetw) {
	//					minBetw=betw2[best];
	//				linkToDelete=i;
	//				}	
	//				adjacency[best][i]=1;
	//				adjacency[i][best]=1;
					
	//	}
	//}

	//adjacency[best][linkToDelete]=0;
	//adjacency[linkToDelete][best]=0;
	
	bool oneFound=false;
	for (int i=0;i<array_size;i++) {
		if ((adjacency[best][i]==1)&&(best!=i)&&(oneFound==false)){
					adjacency[best][i]=0;
					adjacency[i][best]=0;	
					oneFound=true;
		}
	}

	
	//cout<<"The highest value of betweenness is "<<betw[best]<< " (node "<<best+1<<")\n";
}

//Splits the network only if the modularity of the division is higher of the given threshold
//Returns the new q if the splitting has been made, -1 otherwise
double splitNetwork_Threshold (int** adjacency, double *betw, int array_size, double modThreshold) {

	int numberOfGroups_before=0;
	int numberOfGroups_after=0;
	double numberOfLinks_before=0;

	int **groups_before=initialise_int_array(array_size);
	int **groups_after=initialise_int_array(array_size);
	
	double **totalResults=initialise_double_array(array_size);
	int **adjacency_before=initialise_int_array(array_size);	
	
	//adjacency_before is the connectivity matrix before the splitting procedure
	for (int i=0;i<array_size;i++)
		for (int j=0;j<array_size;j++) {
			adjacency_before[i][j]=adjacency[i][j];
		}
	int numberOfMembers_before[array_size];
	int numberOfMembers_after[array_size];

	for (int i=0;i<array_size; i++) {
		numberOfMembers_before[i]=0;
		numberOfMembers_after[i]=0;
	}
	 
	numberOfGroups_before=getGroups(adjacency, groups_before,numberOfMembers_before,array_size);
	
	splitNetwork(adjacency, betw, array_size); 
	
	//adjacency is not renamed but it must be considered adjacency_after
	
	numberOfGroups_after=getGroups(adjacency, groups_after ,numberOfMembers_after,array_size);
	
	double **modularity=initialise_double_array(array_size);
	double **modularity_num=initialise_double_array(array_size);
	
	for (int i=0;i<numberOfGroups_after;i++) 
		for (int j=0;j<numberOfGroups_after;j++) {
			modularity_num[i][j]=0;		
		}
		
	for (int i=0;i<array_size;i++)
		for (int j=0;j<array_size;j++) {
	     if (adjacency_before[i][j]==1)
			  numberOfLinks_before=numberOfLinks_before+1;	
	}


	for (int i=0;i<array_size;i++) 
		for (int j=0;j<array_size;j++) {
			if (adjacency_before[i][j]==1)  {
				for (int g1=0;g1<numberOfGroups_after;g1++) {
					if (isInGroup(i+1,groups_after[g1],numberOfMembers_after[g1])) {
						for (int g2=0;g2<numberOfGroups_after;g2++)
							if (isInGroup(j+1,groups_after[g2],numberOfMembers_after[g2])) {
							
								modularity_num[g1][g2]=modularity_num[g1][g2]+1;
							}
					}
				}
			}
	}		
	
	for (int i=0;i<numberOfGroups_after;i++) 
		for (int j=0;j<numberOfGroups_after;j++) {
			modularity[i][j]=modularity_num[i][j]/numberOfLinks_before;		
		}
   
   //cout<<"The total number of links is:"<<numberOfLinks_before<<"\n";
   //cout<<"The e matrix is the following:\n";
   //print_array_double(modularity_num,numberOfGroups_after);
   
   //calculation of the modularity Q
   //Q= Tr e - ||e^2||
	//where ||x|| indicates the sum of all elements of x  
     
   //calculation of the trace of the matrix e
	double total1=0;   
	for (int i=0;i<numberOfGroups_after;i++)
   	total1=total1+modularity[i][i];
  
	for (int i=0;i<numberOfGroups_after;i++) {
		for (int j=0;j<numberOfGroups_after;j++) {
				double partial=0;
				for (int t=0;t<numberOfGroups_after;t++)
					totalResults[i][j]=partial+modularity[i][t]*modularity[t][j];
		}
	}	  
  
  
	//calculation of ||e^2|| 
	double total2=0;
   for (int i=0;i<numberOfGroups_after;i++)
   	for (int j=0;j<numberOfGroups_after;j++)
   		total2=total2+totalResults[i][j];
   		
   //sum of the two parts Tr e and ||e^2||	
 
   double q=total1-total2;
   
   //cout<<"The scalar value of q is "<<q<<"\n";
   
   
   if (numberOfGroups_after==array_size) {
		//cout<<"It is not possible to split the networks, since the number of groups is equal to the number of hosts.\n";
		return -1;	
	}	
	
	if (q>modThreshold) {
		
		//cout<<"Successful splitting!\n";
		return q;
	
	}

	else {
		for (int i=0;i<array_size;i++)
			for (int j=0;j<array_size;j++) {
				adjacency[i][j]=adjacency_before[i][j];
		}
		return -1;
	}
	
	for (int i=0;i<array_size;i++) {
		delete groups_before[i];
		delete groups_after[i];
		delete totalResults[i];
		delete adjacency_before[i];
		delete modularity[i];
		delete modularity_num[i];
	}
	
	delete groups_before;
	delete groups_after;
	delete totalResults;
	delete adjacency_before;
	delete modularity;
	delete modularity_num;
	
}



//Given the groups of nodes and two nodes "node1" and "node2" returns true if the hosts
//are in the same group, false otherwise
bool areInTheSameGroup (int node1, int node2, int** groups, int numberOfGroups, int* numberOfMembers) {
	
	bool result=false;
	for (int k=0;k<numberOfGroups;k++) {
		
		if (isInGroup(node1,groups[k],numberOfMembers[k]))
			if (isInGroup(node2,groups[k],numberOfMembers[k])) {
				result=true;
				//break;
			}
		}
	
	return result;
}


//print the composition of the groups	
void printGroups(int numberOfGroups, int**groups, int*numberOfMembers, int array_size) {
	
	for (int i=0;i<numberOfGroups;i++) {
		cout<<"\n";
		cout<<"The members of group "<<i+1<<" are: ";
		for (int j=0;j<numberOfMembers[i];j++)
			cout<<groups[i][j]<<" ";
		cout<<"\n";
	}	
}


/*
int main(int argc, char* argv[])
{	int i, j;
	int k;
	double p;
	int    **adjacency 			    	= initialise_int_array(array_sz);
	double **weight					    = initialise_double_array(array_sz);
	double **shortest_path_hops			= initialise_double_array(array_sz);
	double **shortest_path_ideal_hops	= initialise_double_array(array_sz);	
	double **shortest_path_weights		= initialise_double_array(array_sz);
	double **shortest_path_ideal_weights= initialise_double_array(array_sz);
	
	int **groups                        =initialise_int_array(array_sz);
	
	const char *filename = "results.csv";
	fout.open(filename, ios::out | ios::app);
	srand(32);
	
	//experimental area
	//int size=9;
	double threshold=0.05;
	//double **weight=initialise_double_array(size);
	char* filename1=argv[1];
	
	int size=load_double_array(weight,filename1);
	cout<<"\nWeight matrix \n";
	print_array_double(weight,size);
	
	generate_adjacency(weight,adjacency,threshold,size);
	cout<<"\nAdjacency matrix \n";
	print_array_int(adjacency,size);	
	
	calculate_shortest_paths(shortest_path_hops,shortest_path_ideal_hops,adjacency,weight,size);
	cout<<"\nShortest path matrix - number of hops\n";
	print_array_double(shortest_path_hops,size);
	
	cout<<"\nShortest path (ideal) matrix - number of hops\n";
	print_array_double(shortest_path_ideal_hops,size);
		
	calculate_shortest_paths(shortest_path_weights,shortest_path_ideal_weights,adjacency,weight,size);
	cout<<"\nShortest path matrix - weights\n";
	print_array_double(shortest_path_weights,size);
	
	cout<<"\nShortest path (ideal) matrix - weights \n";
	print_array_double(shortest_path_ideal_weights,size);
	
	//calculate betweenness
	double betw[size];
	
	//calculate_betweenness(betw,adjacency,size);
	
	//print betweenness
	//cout<<"\n";
	//for (int i=0; i<size; i++) 
	//	cout<<"Betwenneess of node "<<i+1<<" is "<<betw[i]<<"\n";
	
	//number of the members of the possible groups
	int numberOfMembers[size];
	int numberOfGroups=0;
	//for (int i=0;i<size; i++)
	//	numberOfMembers[i]=0;
	
	//int numberOfGroups=getGroups (adjacency, groups ,numberOfMembers,size);
	//printGroups(numberOfGroups, groups, numberOfMembers, size);

	cout<<"\nAdjacency matrix BEFORE splitting\n";
		print_array_int(adjacency,size);	


	for (int i=0;i<size; i++)
		numberOfMembers[i]=0;
	//splitNetwork(adjacency, betw, size);
 	 	
 	//modularity threshold
 	double modth=0;
 
	bool splitted;
	
	do {
	   splitted=false;
		calculate_betweenness(betw,adjacency,size);
		for (int i=0;i<size; i++)
			numberOfMembers[i]=0;
		   numberOfGroups=getGroups(adjacency, groups ,numberOfMembers,size);
	      printGroups(numberOfGroups, groups, numberOfMembers, size);	
		   modth=splitNetwork_Threshold(adjacency, betw, size, modth);
	}
	while (modth>0);

	cout<<"\nAdjacency matrix AFTER splitting\n";
		print_array_int(adjacency,size);	

	// Set initial lattice parameter
	//k = 10;
	//p=0.2;
	//initialise_adjacency_array(adjacency, array_sz, k, p);
	//initialise_weight_array(weight, array_sz);
	//rewire(adjacency, weight, array_sz, k, p);
	//print_int_array(adjacency, array_sz);
	//print_double_array(weight, array_sz);
	//calculate_shortest_paths(shortest_path, shortest_path_ideal, adjacency, weight, array_sz);

	//print_double_array(shortest_path, array_sz);
	//print_double_array(shortest_path_ideal, array_sz);
	// Work out E_glob
	//
	//double E_glob = calculate_E_glob(shortest_path, shortest_path_ideal, array_sz);
	//double E_loc  = calculate_E_loc(adjacency, weight, array_sz);
	//double cost   = calculate_cost(adjacency, weight, array_sz);
	//fout << p << ", " << k << ", " << E_glob << ", " << E_loc << ", " << cost << "\n"; 
	//fout.flush();
	//cout << p << ", " << k << ", " << E_glob << ", " << E_loc << ", " << cost << "\n"; 
	//cout.flush();


	for (i = 0; i < array_sz; i++)
	{	delete adjacency[i];
		delete weight[i];
		delete shortest_path_hops[i];
		delete shortest_path_ideal_hops[i];
		delete shortest_path_weights[i];
		delete shortest_path_ideal_weights[i];
		delete groups[i];
	}

	delete adjacency;
	delete weight;
	delete shortest_path_hops;
	delete shortest_path_ideal_hops;
	delete shortest_path_weights;
	delete shortest_path_ideal_weights;
	delete groups;
	
	fout.close();

	return 0;
}
*/
