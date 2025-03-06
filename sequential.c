// Sequential version of CM30225 shared memory coursework
// Use the following to compile: gcc -pthread -Wall -Wextra -Wconversion -o sequential sequential.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

// Values MYARRAYSIZE and PRECISION can be changed freely
#define MYARRAYSIZE 100
#define PRECISION 0.01

double** SequentialAlg(double relaxationArray[MYARRAYSIZE][MYARRAYSIZE]);
double findGreatestDifference(double relaxationArray[MYARRAYSIZE][MYARRAYSIZE],
							  double newRelaxationArray[MYARRAYSIZE][MYARRAYSIZE]);
void free2DArray(double **my2DArray);

int main(void) {
	// Variable used for calculating runtime
	clock_t time;
	time = clock();

	int numIterations = 0;

	double relaxationArray[MYARRAYSIZE][MYARRAYSIZE];

	// Declare array values
	for (int y = 0; y < MYARRAYSIZE; y++) {
		for (int x = 0; x < MYARRAYSIZE; x++) {
			relaxationArray[y][x] = rand() % 100;
		}
	}

	// Commented for testing speed
	// printf("Initial Array: \n");
	// for (int y = 0; y < MYARRAYSIZE; y++) {
	// 	for (int x = 0; x < MYARRAYSIZE; x++) {
	// 		printf("%lf ", relaxationArray[y][x]);
	// 	}
	// 	printf("\n");
	// }

	// printf("\n");
	// printf("\n");


	// Set flag for iterating through algorithm until precision is reached
	bool isFinished = false;

	while (isFinished == false) {
		numIterations++;

		double** newRelaxationArray = SequentialAlg(relaxationArray);
		double nonStackRelaxationArray[MYARRAYSIZE][MYARRAYSIZE];

		// Copy results off of the stack
		for (int y = 0; y < MYARRAYSIZE; y++) {
			for (int x = 0; x < MYARRAYSIZE; x++) {
				nonStackRelaxationArray[y][x] = newRelaxationArray[y][x];
			}
		}

		free2DArray(newRelaxationArray);

		// Check if last iteration met the precision criteria
		double myGreatestChange = findGreatestDifference(relaxationArray, nonStackRelaxationArray);
		// Commented for testing speed
		// printf("Change between iterations: %lf \n", myGreatestChange);
		if (myGreatestChange < PRECISION) {
			isFinished = true;
		}

		// Copy for next iteration
		for (int y = 0; y < MYARRAYSIZE; y++) {
			for (int x = 0; x < MYARRAYSIZE; x++) {
				relaxationArray[y][x] = nonStackRelaxationArray[y][x];
			}
		}
	}

	// Commented for testing speed
	// printf("\n Final Array: \n");
	// for (int y = 0; y < MYARRAYSIZE; y++) {
	// 	for (int x = 0; x < MYARRAYSIZE; x++) {
	// 		printf("%lf ", relaxationArray[y][x]);
	// 	}
	// 	printf("\n");
	// }

	// printf("\n");
	// printf("\n");

	// Get time taken in seconds
	time = clock() - time;
	double timeTaken = ((double)time)/CLOCKS_PER_SEC;

	printf("Program took %f seconds and %d iterations to execute \n", timeTaken, numIterations);

	return 0;
}

// Carries out the relaxation algorithm on a split 2D array sequentially
// Takes in a 2d array of fixed size (array[double][double])
// Returns nothing
// ALLOCATES MEMORY
double** SequentialAlg(double relaxationArray[MYARRAYSIZE][MYARRAYSIZE]) {
	// Creates a new relaxation array for returning, and allocates memory to it
	double** newRelaxationArray = (double**)malloc(MYARRAYSIZE * sizeof(double*));
    for (int y = 0; y < MYARRAYSIZE; y++) {
        newRelaxationArray[y] = (double*)malloc(MYARRAYSIZE * sizeof(double));
    }

	for (int y = 0; y < MYARRAYSIZE; y++) {
		for (int x = 0; x < MYARRAYSIZE; x++) {
			newRelaxationArray[y][x] = relaxationArray[y][x];
		}
	}

	// Average non-border elements
	for (int y = 1; y < MYARRAYSIZE - 1; y++) {
		for (int x = 1; x < MYARRAYSIZE - 1; x++) {
			newRelaxationArray[y][x] = (relaxationArray[y][x + 1] + 
			                            relaxationArray[y][x - 1] + 
			                            relaxationArray[y + 1][x] + 
			                            relaxationArray[y - 1][x]) / 4;

		}
	}

	return newRelaxationArray;
}

// Finds the greatest change between two 2D arrays
// Takes in the original and new relaxationArray (double[MYARRAYSIZE][MYARRAYSIZE])
// Returns the greatest change between the arrays as a double
double findGreatestDifference(double relaxationArray[MYARRAYSIZE][MYARRAYSIZE],
							  double newRelaxationArray[MYARRAYSIZE][MYARRAYSIZE]) {
	double greatestChange = 0.0;
	double currentChange;

	for (int y = 0; y < MYARRAYSIZE; y++) {
		for (int x = 0; x < MYARRAYSIZE; x++) {
			currentChange = relaxationArray[y][x] - newRelaxationArray[y][x];

			if (currentChange > greatestChange) {
				greatestChange = currentChange;
			}
		}
	}

	return greatestChange;
}

// Frees allocated memory
void free2DArray(double **my2DArray) {
	for (int y = 0; y < MYARRAYSIZE; y++) {
		free(my2DArray[y]);
	}
	free(my2DArray);
}