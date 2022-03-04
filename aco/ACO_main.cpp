#include "ACO.h"

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <limits>
#include <climits>
#include <fstream>

using namespace std;

ACO::ACO (int nAnts, int nCities, 
		double alpha, double beta, double q, double ro, double taumax,
		int initCity) {
	NUMBEROFANTS 	= nAnts;
	NUMBEROFCITIES 	= nCities;
	ALPHA 			= alpha;
	BETA 			= beta;
	Q 				= q;
	RO 				= ro;
	TAUMAX 			= taumax;
	INITIALCITY		= initCity;

	randoms = new Randoms (21);	
}
ACO::~ACO () {
	for(int i=0; i<NUMBEROFCITIES; i++) {
		delete [] GRAPH[i];
		delete [] CITIES[i];
		delete [] PHEROMONES[i];
		delete [] DELTAPHEROMONES[i];
		if(i < NUMBEROFCITIES - 1) {
			delete [] PROBS[i];
		}
	}
	delete [] GRAPH;
	delete [] CITIES;
	delete [] PHEROMONES;
	delete [] DELTAPHEROMONES;
	delete [] PROBS;
}

void ACO::init () {
	GRAPH 			= new int*[NUMBEROFCITIES];
	CITIES 			= new double*[NUMBEROFCITIES];
	PHEROMONES 		= new double*[NUMBEROFCITIES];
	DELTAPHEROMONES = new double*[NUMBEROFCITIES];
	PROBS 			= new double*[NUMBEROFCITIES-1];
	for(int i=0; i<NUMBEROFCITIES; i++) {
		GRAPH[i] 			= new int[NUMBEROFCITIES];
		CITIES[i] 			= new double[2];
		PHEROMONES[i] 		= new double[NUMBEROFCITIES];
		DELTAPHEROMONES[i] 	= new double[NUMBEROFCITIES];
		PROBS[i] 			= new double[2];
		for (int j=0; j<2; j++) {
			CITIES[i][j] = -1.0;
			PROBS[i][j]  = -1.0;
		}
		for (int j=0; j<NUMBEROFCITIES; j++) {
			GRAPH[i][j] 			= 0;
			PHEROMONES[i][j] 		= 0.0;
			DELTAPHEROMONES[i][j] 	= 0.0;
		}
	}	

	ROUTES = new int*[NUMBEROFANTS];
	for (int i=0; i<NUMBEROFANTS; i++) {
		ROUTES[i] = new int[NUMBEROFCITIES];
		for (int j=0; j<NUMBEROFCITIES; j++) {
			ROUTES[i][j] = -1;
		}
	}
	
	BESTLENGTH = (double) INT_MAX;
	BESTROUTE  = new int[NUMBEROFCITIES];
	for (int i=0; i<NUMBEROFCITIES; i++) {
		BESTROUTE[i] = -1;	
	}
}


void ACO::connectCITIES (int cityi, int cityj) {
	GRAPH[cityi][cityj] = 1;
	PHEROMONES[cityi][cityj] = randoms -> Uniforme() * TAUMAX;
	GRAPH[cityj][cityi] = 1;
	PHEROMONES[cityj][cityi] = PHEROMONES[cityi][cityj];
}
void ACO::setCITYPOSITION (int city, double x, double y) {
	CITIES[city][0] = x;
	CITIES[city][1] = y;
}
void ACO::printPHEROMONES () {	
	cout << " PHEROMONES: " << endl;
	cout << "  | ";
	for (int i=0; i<NUMBEROFCITIES; i++) {
		printf("%5d   ", i);
	}
	cout << endl << "- | ";
	for (int i=0; i<NUMBEROFCITIES; i++) {
		cout << "--------";
	}
	cout << endl;
	for (int i=0; i<NUMBEROFCITIES; i++) {
		cout << i << " | ";
		for (int j=0; j<NUMBEROFCITIES; j++) {
			if (i == j) {
				printf ("%5s   ", "x");
				continue;
			}
			if (exists(i, j)) {
				printf ("%7.3f ", PHEROMONES[i][j]);
			}
			else {
				if(PHEROMONES[i][j] == 0.0) {
					printf ("%5.0f   ", PHEROMONES[i][j]);
				}
				else {
					printf ("%7.3f ", PHEROMONES[i][j]);
				}
			}
		}
		cout << endl;
	}
	cout << endl;
}


int ACO::distance (int cityi, int cityj) {
	return (int) 
		sqrt (pow (CITIES[cityi][0] - CITIES[cityj][0], 2) + 
 			  pow (CITIES[cityi][1] - CITIES[cityj][1], 2) + 0.5);
}
bool ACO::exists (int cityi, int cityc) {
	return (GRAPH[cityi][cityc] == 1);
}
bool ACO::vizited (int antk, int c) {
	for (int l=0; l<NUMBEROFCITIES; l++) {
		if (ROUTES[antk][l] == -1) {
			break;
		}
		if (ROUTES[antk][l] == c) {
			return true;
		}
	}
	return false;
}
double ACO::PHI (int cityi, int cityj, int antk) {
	double ETAij = (double) pow (1.0 / distance (cityi, cityj), BETA);
	double TAUij = (double) pow (PHEROMONES[cityi][cityj],   ALPHA);

	double sum = 0.0;
	for (int c=0; c<NUMBEROFCITIES; c++) {
		if (exists(cityi, c)) {
			if (!vizited(antk, c)) {
				double ETA = (double) pow (1.0 / distance (cityi, c), BETA);
				double TAU = (double) pow (PHEROMONES[cityi][c],   ALPHA);
				sum += ETA * TAU;
			}	
		}	
	}
	return (ETAij * TAUij) / sum;
}

double ACO::length (int antk) {
	double sum = 0.0;
	for (int j=0; j<NUMBEROFCITIES-1; j++) {
		sum += distance (ROUTES[antk][j], ROUTES[antk][j+1]);	
	}
	return sum;
}

int ACO::city () {
	double xi = randoms -> Uniforme();
	int i = 0;
	double sum = PROBS[i][0];
	while (sum < xi) {
		i++;
		sum += PROBS[i][0];
	}
	return (int) PROBS[i][1];
}

void ACO::route (int antk) {
	ROUTES[antk][0] = INITIALCITY;
	for (int i=0; i<NUMBEROFCITIES-1; i++) {		
		int cityi = ROUTES[antk][i];
		int count = 0;
		for (int c=0; c<NUMBEROFCITIES; c++) {
			if (cityi == c) {
				continue;	
			}
			if (exists (cityi, c)) {
				if (!vizited (antk, c)) {
					PROBS[count][0] = PHI (cityi, c, antk);
					PROBS[count][1] = (double) c;
					count++;
				}

			}
		}
		
		// deadlock
		if (0 == count) {
			return;
		}
		
		ROUTES[antk][i+1] = city();
	}
}
int ACO::valid (int antk, int iteration) {
	for(int i=0; i<NUMBEROFCITIES-1; i++) {
		int cityi = ROUTES[antk][i];
		int cityj = ROUTES[antk][i+1];
		if (cityi < 0 || cityj < 0) {
			return -1;	
		}
		if (!exists(cityi, cityj)) {
			return -2;	
		}
		for (int j=0; j<i-1; j++) {
			if (ROUTES[antk][i] == ROUTES[antk][j]) {
				return -3;
			}	
		}
	}
	
	if (!exists (INITIALCITY, ROUTES[antk][NUMBEROFCITIES-1])) {
		return -4;
	}
	
	return 0;
}

void ACO::printGRAPH () {
	cout << " GRAPH: " << endl;
	cout << "  | ";
	for( int i=0; i<NUMBEROFCITIES; i++) {
		cout << i << " ";
	}
	cout << endl << "- | ";
	for (int i=0; i<NUMBEROFCITIES; i++) {
		cout << "- ";
	}
	cout << endl;
	int count = 0;
	for (int i=0; i<NUMBEROFCITIES; i++) {
		cout << i << " | ";
		for (int j=0; j<NUMBEROFCITIES; j++) {
			if(i == j) {
				cout << "x ";	
			}
			else {
				cout << GRAPH[i][j] << " ";	
			}
			if (GRAPH[i][j] == 1) {
				count++;	
			}
		}
		cout << endl;
	}
	cout << endl;
//	cout << "Number of connections: " << count << endl << endl;
}
void ACO::printRESULTS () {
	BESTLENGTH += distance (BESTROUTE[NUMBEROFCITIES-1], INITIALCITY);
	cout << " BEST ROUTE:" << endl;
	for (int i=0; i<NUMBEROFCITIES; i++) {
		cout << BESTROUTE[i] << " ";
	}
	cout << endl << "length: " << BESTLENGTH << endl;
	
	cout << endl << " IDEAL ROUTE:" << endl;
	// cout << "0 7 6 2 4 5 1 3" << endl;
	// cout << "length: 127.509" << endl;
}

void ACO::updatePHEROMONES () {
	for (int k=0; k<NUMBEROFANTS; k++) {
		double rlength = length(k);
		for (int r=0; r<NUMBEROFCITIES-1; r++) {
			int cityi = ROUTES[k][r];
			int cityj = ROUTES[k][r+1];
			DELTAPHEROMONES[cityi][cityj] += Q / rlength;
			DELTAPHEROMONES[cityj][cityi] += Q / rlength;
		}
	}
	for (int i=0; i<NUMBEROFCITIES; i++) {
		for (int j=0; j<NUMBEROFCITIES; j++) {
			PHEROMONES[i][j] = (1 - RO) * PHEROMONES[i][j] + DELTAPHEROMONES[i][j];
			DELTAPHEROMONES[i][j] = 0.0;
		}	
	}
}


void ACO::optimize (int ITERATIONS) {
	for (int iterations=1; iterations<=ITERATIONS; iterations++) {
		cout << flush;
//		cout << "ITERATION " << iterations << " HAS STARTED!" << endl << endl;

		for (int k=0; k<NUMBEROFANTS; k++) {
//			cout << " : ant " << k << " has been released!" << endl;
			while (0 != valid(k, iterations)) {
//				cout << "  :: releasing ant " << k << " again!" << endl;
				for (int i=0; i<NUMBEROFCITIES; i++) {
					ROUTES[k][i] = -1;	
				}
				route(k);
			}
			
			// for (int i=0; i<NUMBEROFCITIES; i++) {
			// 	cout << ROUTES[k][i] << " ";	
			// }
			// cout << endl;
			
//			cout << "  :: route done" << endl;
			double rlength = length(k);

			if (rlength < BESTLENGTH) {
				BESTLENGTH = rlength;
				for (int i=0; i<NUMBEROFCITIES; i++) {
					BESTROUTE[i] = ROUTES[k][i];
				}
			}
//			cout << " : ant " << k << " has ended!" << endl;				
		}		

//		cout << endl << "updating PHEROMONES . . .";
		updatePHEROMONES ();
//		cout << " done!" << endl << endl;
//		printPHEROMONES ();
		
		for (int i=0; i<NUMBEROFANTS; i++) {
			for (int j=0; j<NUMBEROFCITIES; j++) {
				ROUTES[i][j] = -1;
			}
		}

//		cout << endl << "ITERATION " << iterations << " HAS ENDED!" << endl << endl;
	}
}



int ITERATIONS;

int NUMBEROFANTS;
int NUMBEROFCITIES;

// if (ALPHA == 0) { stochastic search & sub-optimal route }
#define ALPHA			(double) 0.5
// if (BETA  == 0) { sub-optimal route }
#define BETA			(double) 0.8
// Estimation of the suspected best route.
#define Q				(double) 80
// Pheromones evaporation. 
#define RO				(double) 0.2
// Maximum pheromone random number.
#define TAUMAX			(int) 2

#define INITIALCITY		(int) 0
#define OUT "aco_out_15.txt"
#define IN  "input15.txt"


int main() {

	ifstream infile;
	infile.open(IN);
	infile >> NUMBEROFCITIES;

	NUMBEROFANTS = ITERATIONS = NUMBEROFCITIES * 8 / 10;
	ACO *ANTS = new ACO (NUMBEROFANTS, NUMBEROFCITIES, 
			 			ALPHA, BETA, Q, RO, TAUMAX,
			 			INITIALCITY);
	ITERATIONS += 50 * NUMBEROFCITIES;
	ANTS -> init();

	for(int i = 0; i < NUMBEROFCITIES; i++){
		for(int j = 0; j < i; j++) ANTS -> connectCITIES(i, j);
	}
	/* ANTS -> connectCITIES(8, 2);
	ANTS -> connectCITIES(8, 6);
	ANTS -> connectCITIES(8, 7); */

	ANTS -> setCITYPOSITION (0,  1,  1);
	for(int i = 0; i < NUMBEROFCITIES; i++){
		int x, y;
		infile >> x >> y;
		ANTS -> setCITYPOSITION (i, x, y);
	}
	infile.close();
	// ANTS -> setCITYPOSITION(8, 26, 20);
	ofstream outfile;
	outfile.open(OUT);
	clock_t st = clock();
//	ANTS -> printGRAPH ();

//	ANTS -> printPHEROMONES ();

	ANTS -> optimize (ITERATIONS);

	ANTS -> printRESULTS ();

	clock_t end =  clock();
	cout << "\nTime of the algo: " <<  ((double) (end - st) / CLOCKS_PER_SEC) << '\n';
	outfile << "Time of the algo: " <<  ((double) (end - st) / CLOCKS_PER_SEC) << '\n';
	outfile << "length: " << ANTS -> BESTLENGTH << '\n';
	// for(int i = 0; i < NUMBEROFCITIES; i++){
	// 	for(int j = 0; j < NUMBEROFCITIES; j++) cout << ANTS -> distance(i, j) << ' ';
	// 	cout << '\n';
	// }
	return 0;
}
