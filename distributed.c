// Multicore version of CM30225 distributed memory coursework
// Use the following to compile: mpicc -Wall -Wconversion -Wextra -o distributed distributed.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


void sequentialAlg(int segmentRows, int largestSegmentRows, int segmentColumns, 
				   double relaxationArray[segmentRows][segmentColumns]);
void findGreatestDifference(double greatestChange[1], int segmentRows, int segmentColumns,
							  double relaxationArray[segmentRows][segmentColumns],
							  double newRelaxationArray[segmentRows][segmentColumns]);
void joinArrays(int numNodes, int firstSegmentRows, int segmentRows, int finalSegmentRows,
                int largestSegmentRows, int mainArraySize,
				double relaxationArray[largestSegmentRows][mainArraySize],
                double allRelaxationArrays[numNodes][largestSegmentRows][mainArraySize]);


// Following arguments are to be passed in as parameters:
// mainArraySize - Length and width of square array
// precision - Precision to be relaxed to until no changes of this size are made 
// (e.g. 3 gives precision 0.001)
int main(int argc, char **argv) {
	int mainArraySize = atoi(argv[1]);
	int precisionPower = atoi(argv[2]);
	double precision = pow(0.1, precisionPower);


	int rc, myRank, numNodes;
	rc = MPI_Init(&argc, &argv);

	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numNodes);

	// Calculate runtime and num iterations using master node
	clock_t time;
	if(myRank == 0) {
		time = clock();
	}
	int numIterations = 0;

	// Must account for the fact that the root node doesn't do normal processing 
	int segmentRows = (mainArraySize / (numNodes - 1)) + 2;
	int firstSegmentRows = segmentRows - 1;
	int finalSegmentRows = firstSegmentRows + (mainArraySize % (numNodes - 1));
	int largestSegmentRows;
	if (finalSegmentRows > segmentRows) {
		largestSegmentRows = finalSegmentRows;
	} else {
		largestSegmentRows = segmentRows;
	}

	// Segments are smaller portions of the main array to be processed by a node
	double relaxationArray[mainArraySize][mainArraySize];
	double segmentArray[largestSegmentRows][mainArraySize];
	double oldSegmentArray[largestSegmentRows][mainArraySize];
	// 3D array is declared for gathering arrays from other nodes
	double allSegmentArrays[numNodes][largestSegmentRows][mainArraySize];
	// 3D array for sending first and last array of each segment
	double allSndSndLastArrays[numNodes][2][mainArraySize];
	double allFirstLastArrays[numNodes][2][mainArraySize];
	// 2D arrays for containing these values for each segment
	double sndSndLastArrays[2][mainArraySize];
	double firstLastArrays[2][mainArraySize];
	
	// Used for calculating precision, needs to be an array for use with MPI
	double greatestChange[1] = {0.11};
	int finished[1] = {0};

	// Process with rank 0 is master, so has very different actions to take
	if (myRank == 0) {
		// Start by broadcasting the initial relaxation array to every slave node
		for (int y = 0; y < mainArraySize; y++) {
			for (int x = 0; x < mainArraySize; x++) {
				relaxationArray[y][x] = rand() % 100;
			}
		}

		MPI_Bcast(relaxationArray, mainArraySize * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// Iterate until defined precision is reached
		while (finished[0] == 0) {
			// Gather the maximum change across all segments for deciding when to finish
			greatestChange[0] = 0;
			MPI_Reduce(greatestChange, greatestChange, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			if (greatestChange[0] < precision) {
				finished[0] = 1;
			}

			// Broadcast if finished to other nodes
			MPI_Bcast(finished, 1, MPI_INT, 0, MPI_COMM_WORLD);

			if (finished[0] == 0) {
				numIterations++;

				// Gather the second and second to last arrays from each segment
				MPI_Gather(sndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, 
				           allSndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Arrange these arrays to find the border arrays for other segments that must be replaced
				for (int x = 0; x < mainArraySize; x++) {
					allFirstLastArrays[1][1][x] = allSndSndLastArrays[2][0][x];
					allFirstLastArrays[numNodes - 1][0][x] = allSndSndLastArrays[numNodes - 2][1][x];
				}

				for (int z = 2; z < numNodes - 1; z++) {
					for (int x = 0; x < mainArraySize; x++) {
						allFirstLastArrays[z][0][x] = allSndSndLastArrays[z - 1][1][x];
						allFirstLastArrays[z][1][x] = allSndSndLastArrays[z + 1][0][x];
					}
				}

				// Send the updated first and last array of each segment out
				MPI_Scatter(allFirstLastArrays, 2 * mainArraySize, MPI_DOUBLE,
							firstLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// Otherwise, collect all arrays and rejoin to produce final array
			} else {
				MPI_Gather(segmentArray, largestSegmentRows * mainArraySize, MPI_DOUBLE, 
				           allSegmentArrays, largestSegmentRows * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				joinArrays(numNodes, firstSegmentRows, segmentRows, finalSegmentRows, largestSegmentRows,
				           mainArraySize, relaxationArray, allSegmentArrays);
			}
		}
	
	// Once a solution has been found, get runtime and terminate all processes
	// Get time taken in seconds
	time = clock() - time;
	double timeTaken = ((double)time)/CLOCKS_PER_SEC;

	printf("Program took %f seconds and %d iterations to execute \n", timeTaken, numIterations);

	MPI_Abort(MPI_COMM_WORLD, rc);

	
	// First and last slave process have slightly different behaviour due to border values
	// Code reuse here, should refactor
	} else if (myRank == 1) {
		// Start by receiving the initial array, and copying the segment over
		MPI_Bcast(relaxationArray, mainArraySize * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int y = 0; y < firstSegmentRows; y++) {
			for (int x = 0; x < mainArraySize; x++) {
				segmentArray[y][x] = relaxationArray[y][x];
				oldSegmentArray[y][x] = oldSegmentArray[y][x];
			}
		}

		// Iterate until defined precision is reached
		while (finished[0] == 0) {
			// Send the maximum change for the segment, for deciding when to finish
			MPI_Reduce(greatestChange, greatestChange, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			// Receive whether or not the relaxation is finished
			MPI_Bcast(finished, 1, MPI_INT, 0, MPI_COMM_WORLD);

			// If it isn't, continue relaxation
			if (finished[0] == 0) {
				// Perform processing on array
				sequentialAlg(firstSegmentRows, largestSegmentRows, mainArraySize, segmentArray);

				// Send the second to last array from the segment (second array here is unnecessary)
				for (int x = 0; x < mainArraySize; x++) {
					sndSndLastArrays[1][x] = segmentArray[firstSegmentRows - 2][x];
				}

				MPI_Gather(sndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, allSndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Receive the new first and last array, from another node
				MPI_Scatter(allFirstLastArrays, 2 * mainArraySize, MPI_DOUBLE,
							firstLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				for (int x = 0; x < mainArraySize; x++) {
					segmentArray[firstSegmentRows - 1][x] = firstLastArrays[1][x];
				}

				// Caclulate greatest change
				findGreatestDifference(greatestChange, firstSegmentRows, mainArraySize, oldSegmentArray, segmentArray);
				for (int y = 0; y < largestSegmentRows; y++) {
					for (int x = 0; x < mainArraySize; x++) {
						oldSegmentArray[y][x] = segmentArray[y][x];
					}
				}

			// Otherwise, send current array to root node
			} else {
				MPI_Gather(segmentArray, largestSegmentRows * mainArraySize, MPI_DOUBLE, 
				           allSegmentArrays, largestSegmentRows * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
		}
	} else if (myRank == (numNodes - 1)) {
		// Start by receiving the initial array, and copying the segment over
		MPI_Bcast(relaxationArray, mainArraySize * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int y = 0; y < finalSegmentRows; y++) {
			for (int x = 0; x < mainArraySize; x++) {
				segmentArray[y][x] = relaxationArray[y + (mainArraySize - finalSegmentRows)][x];
				oldSegmentArray[y][x] = segmentArray[y][x];
			}
		}

		// Iterate until defined precision is reached
		while (finished[0] == 0) {
			// Send the maximum change for the segment, for deciding when to finish
			MPI_Reduce(greatestChange, greatestChange, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			// Receive whether or not the relaxation is finished
			MPI_Bcast(finished, 1, MPI_INT, 0, MPI_COMM_WORLD);

			// If it isn't, continue relaxation
			if (finished[0] == 0) {
				// Perform processing on array
				sequentialAlg(finalSegmentRows, largestSegmentRows, mainArraySize, segmentArray);

				// Send the second array from the segment (second to last array here is unnecessary)
				for (int x = 0; x < mainArraySize; x++) {
					sndSndLastArrays[0][x] = segmentArray[1][x];
				}

				MPI_Gather(sndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, allSndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Receive the new first and last array, from another node
				MPI_Scatter(allFirstLastArrays, 2 * mainArraySize, MPI_DOUBLE,
							firstLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				for (int x = 0; x < mainArraySize; x++) {
					segmentArray[0][x] = firstLastArrays[0][x];
				}

				// Caclulate greatest change
				findGreatestDifference(greatestChange, finalSegmentRows, mainArraySize, oldSegmentArray, segmentArray);
				for (int y = 0; y < largestSegmentRows; y++) {
					for (int x = 0; x < mainArraySize; x++) {
						oldSegmentArray[y][x] = segmentArray[y][x];
					}
				}

			// Otherwise, send current array to root node
			} else {
				MPI_Gather(segmentArray, largestSegmentRows * mainArraySize, MPI_DOUBLE, 
				           allSegmentArrays, largestSegmentRows * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
		}
	// And for all other nodes
	} else {
		// Start by receiving the initial array, and copying the segment over
		MPI_Bcast(relaxationArray, mainArraySize * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int y = 0; y < segmentRows; y++) {
			for (int x = 0; x < mainArraySize; x++) {
				segmentArray[y][x] = relaxationArray[y + ((segmentRows - 2) * (myRank - 1)) - 1][x];
				oldSegmentArray[y][x] = segmentArray[y][x];
			}
		}

		// Iterate until defined precision is reached
		while (finished[0] == 0) {
			// Send the maximum change for the segment, for deciding when to finish
			MPI_Reduce(greatestChange, greatestChange, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			// Receive whether or not the relaxation is finished
			MPI_Bcast(finished, 1, MPI_INT, 0, MPI_COMM_WORLD);

			// If it isn't, continue relaxation
			if (finished[0] == 0) {
				// Perform processing on array
				sequentialAlg(segmentRows, largestSegmentRows, mainArraySize, segmentArray);

				// Send the second and second to last
				for (int x = 0; x < mainArraySize; x++) {
					sndSndLastArrays[0][x] = segmentArray[1][x];
					sndSndLastArrays[1][x] = segmentArray[firstSegmentRows - 1][x];
				}
				MPI_Gather(sndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, allSndSndLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Receive the new first and last array, from another node
				MPI_Scatter(allFirstLastArrays, 2 * mainArraySize, MPI_DOUBLE,
							firstLastArrays, 2 * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				for (int x = 0; x < mainArraySize; x++) {
					segmentArray[0][x] = firstLastArrays[0][x];
					segmentArray[segmentRows - 1][x] = firstLastArrays[1][x];
				}

				// Caclulate greatest change
				findGreatestDifference(greatestChange, segmentRows, mainArraySize, oldSegmentArray, segmentArray);
				for (int y = 0; y < largestSegmentRows; y++) {
					for (int x = 0; x < mainArraySize; x++) {
						oldSegmentArray[y][x] = segmentArray[y][x];
					}
				}

			// Otherwise, send current array to root node
			} else {
				MPI_Gather(segmentArray, largestSegmentRows * mainArraySize, MPI_DOUBLE, 
				           allSegmentArrays, largestSegmentRows * mainArraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
		}
	}
	
	
	MPI_Finalize();
	
	return 0;
}


// Carries out the relaxation algorithm on a split 2D array sequentially
// Takes in 2 2d arrays of fixed size (array[double][double]), one for calculation, and another to store results
// Returns nothing
void sequentialAlg(int segmentRows, int largestSegmentRows, int segmentColumns, 
				   double relaxationArray[largestSegmentRows][segmentColumns]) {

	double newRelaxationArray[largestSegmentRows][segmentColumns];

	// Average non-border elements
	for (int y = 1; y < segmentRows - 1; y++) {
		for (int x = 1; x < segmentColumns - 1; x++) {
			newRelaxationArray[y][x] = (relaxationArray[y][x + 1] + 
			                            relaxationArray[y][x - 1] + 
			                            relaxationArray[y + 1][x] + 
			                            relaxationArray[y - 1][x]) / 4;

		}
	}

	for (int y = 1; y < segmentRows - 1; y++) {
		for (int x = 1; x < segmentColumns - 1; x++) {
			relaxationArray[y][x] = newRelaxationArray[y][x];
		}
	}
}


// Finds the greatest change between two 2D arrays
// Takes in the original and new relaxationArray (double[segmentRows][segmentColumns])
// Returns the greatest change between the arrays
void findGreatestDifference(double greatestChange[1], int segmentRows, int segmentColumns,
							  double relaxationArray[segmentRows][segmentColumns],
							  double newRelaxationArray[segmentRows][segmentColumns]) {
	greatestChange[0] = 0.0;
	double currentChange;

	for (int y = 0; y < segmentRows; y++) {
		for (int x = 0; x < segmentColumns; x++) {
			currentChange = relaxationArray[y][x] - newRelaxationArray[y][x];

			if (currentChange > greatestChange[0]) {
				greatestChange[0] = currentChange;
			}
		}
	}
}


// Join all of the arrays gathered from each node to one full array after each iteration
// Takes in the target array and array of gathered arrays
// Returns no input
void joinArrays(int numNodes, int firstSegmentRows, int segmentRows, int finalSegmentRows,
                int largestSegmentRows, int mainArraySize,
				double relaxationArray[largestSegmentRows][mainArraySize],
                double allRelaxationArrays[numNodes][largestSegmentRows][mainArraySize]) {

	// Start at z = 1 as first array is from root node
	// Only copy the lines that were changed and outer boundaries
	for (int y = 0; y < firstSegmentRows - 1; y++) {
		for (int x = 0; x < mainArraySize; x++) {
			relaxationArray[y][x] = allRelaxationArrays[1][y][x];
		}
	}

	for (int z = 2; z < numNodes - 1; z++) {
		for (int y = 0; y < segmentRows - 2; y++) {
			for (int x = 0; x < mainArraySize; x++) {
				relaxationArray[y + ((z - 1) * (segmentRows - 2))][x] = allRelaxationArrays[z][y + 1][x];
			}
		}
	}

	for (int y = 1; y < finalSegmentRows; y++) {
		for (int x = 0; x < mainArraySize; x++) {
			relaxationArray[y + (mainArraySize - (finalSegmentRows))][x] = allRelaxationArrays[numNodes - 1][y][x];
		}
	}

}