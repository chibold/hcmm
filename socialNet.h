// socialNet.h
//
// A library for the simulation of complex social networks
//
// September 2004  - March 2005
// 
// Author Mirco Musolesi
// 

void print_int_array(int **array, int array_size);

void print_double_array(double **array, int array_size);

void print_float_array(float **array, int array_size);

int **initialise_int_array(int array_size);

double **initialise_double_array(int array_size);

void initialise_adjacency_array_WS(int **adjacency, int n, int k, double p);

void initialise_adjacency_array_random(int **adjacency, int array_size);

void initialise_adjacency_array(int **adjacency, int array_size, int k, double p);

void initialise_weight_array_WS(double **weight, int array_size);

void initialise_weight_array_random(double **weight, int array_size);

void initialise_weight_array_ingroups(double ** weight, int array_size, int numberOfGroups, double nRewiring, double threshold, double seed);

void initialise_weight_array_ingroups2(double ** weight, int array_size, int numberOfGroups, double nRewiring, double threshold, double seed);

void rewire_WS(int **adjacency, double **weight, int n, int k, double p);

void rewire_WS2(int **adjacency, double **weight, int n, int k, double p);

void initialise_weight_array(double **weight, int array_size);

void initialise_weight_array2(double ** weight, int array_size);

void rewire(int **adjacency, double **weight, int array_size, int k, double p);

void calculate_shortest_paths(double **shortest_path, double **shortest_path_ideal, int **adjacency, double **weight, int array_size);

double calculate_energy(double **shortest_path, int array_size);

double calculate_energy_local(int node, int **adjacency, double **weight, int array_size);

double calculate_E_glob(double **shortest_path, double **shortest_path_ideal, int array_size);

double calculate_E_loc(int **adjacency, double **weight, int array_size);

double gamma(double l);

double calculate_cost(int **adjacency, double **weight, int array_size);

int load_double_array(double** array,int array_size, const char* filename);

int load_double_array(double** array, const char* filename);

void generate_adjacency (double** weightMat, int** adjacencyMat, double threshold, int array_size);

void print_array_int (int** array, int array_size);

void print_array_double (double** array, int array_size);

void assign_distance_not_working(int current, int previous, double* distance, int d, int target, int* assigned, int**adjacency, int **pred,int array_size);

void assign_distance(double* distance, int d, int* assigned, int**adjacency, int **pred, int* predNum, int array_size);

void calculate_b (double* betw, double* distance, int **pred, int* predNum, int array_size);

void calculate_betweenness(double* result_betw, int** adjacency, int array_size);

int assignToAGroup (int current, int currentGroup, int** adjacency, int** groups, int* numberOfMembers,int array_size,bool* assigned);

int getGroups (int** adjacency, int** groups, int* numberOfMembers, int array_size);

bool isInGroup(int node, int* group, int numberOfMembers);

double splitNetwork (int** adjacency, double *betw,int array_size);

double splitNetwork_Threshold (int** adjacency, double *betw, int array_size, double modThreshold);

bool areInTheSameGroup (int node1, int node2, int** groups, int numberOfGroups, int* numberOfMembers);

void printGroups(int numberOfGroups, int**groups, int*numberOfMembers, int array_size);
