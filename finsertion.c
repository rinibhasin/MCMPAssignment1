
#include <float.h>
#include <stdlib.h>
#include<stdbool.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
// Declare the functions from coordReader.c file
int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);
void farthestInsertion_TSP(double **distances, int numOfCoords, char *outputfile);
double calculateDistance(double x1, double y1, double x2, double y2);
double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix);
/*void printDistanceMatrix(double **distanceMatrix, int numOfCoords) {
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            printf("%.2f\t", distanceMatrix[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
}*/

int main(int argc, char *argv[]) {
    double start, end;
    double cpu_time_used;
    // Check if the correct number of command-line arguments is provided
    if (argc != 3) {
        printf("Usage: %s arg1 arg2\n", argv[0]);
        return 1; // Return an error code
    }
    // Extract the two command-line arguments
    char *coordinatefile = argv[1];
    char *outputfile = argv[2];
    // Print the arguments
//    printf("Argument 1: %s\n", coordinatefile);
//    printf("Argument 2: %s\n", outputfile);
// Read the coordinates and number of coordinates from .coord file
    int numOfCoords = readNumOfCoords(coordinatefile);
    double **coordinates = readCoords(coordinatefile, numOfCoords);
    double **distanceMatrix = (double **) malloc(numOfCoords * sizeof(double *));
    int i;
    for (i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *) malloc(numOfCoords * sizeof(double));
    }
// Start the clock
    start = omp_get_wtime();
// Calculate the Distance matrix
    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);
// Call the farthestInsertionTSP function to find the tour
    farthestInsertion_TSP(distanceMatrix, numOfCoords, outputfile);
    end = omp_get_wtime();
    cpu_time_used = end - start;
    printf("CPU time used: %f seconds\n", cpu_time_used);
    free(distanceMatrix);
}

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}
double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix) {
    int i, j;
    for (i = 0; i < numOfCoords; i++) {
        for (j = 0; j < numOfCoords; j++) {
            double x1 = coordinates[i][0];
            double y1 = coordinates[i][1];
            double x2 = coordinates[j][0];
            double y2 = coordinates[j][1];
            double distance = calculateDistance(x1, y1, x2, y2);
//                printf("The distance calculated for %f:", distance);
//                printf("\n");
            distanceMatrix[i][j] = distance;
//            printf("%f\t", distanceMatrix[i][j]);
        }
        //      printf("\n");
    }
    return distanceMatrix;
}
// Function to solve the TSP using Farthest insertion method
void farthestInsertion_TSP(double **distances, int numOfCoords, char *outputfile) {
    double minimumDistance = -1;
    int currentSize = 1;
    int *tour = malloc(numOfCoords * sizeof(int));
    int *visited_nodes = malloc(numOfCoords + 2 * sizeof(int));
    int farthest = 0;
    int noOfThreads = omp_get_max_threads();
    int k;
    for (k = 0; k < numOfCoords + 2; k++) {
        visited_nodes[k] = 0;
    }

    /*
        * Step 1 - Start off with a vertex V0
        */
    tour[0] = 0;
    visited_nodes[0] = true;
    int i = 0;
    /*
     * Step 2 -  Find a vertex Vi such that dist(V0, Vi ) is maximal, and create a partial tour (V0, Vi, V0)
     */
    for (i = 1; i < numOfCoords; i++) {
        if (distances[0][i] > minimumDistance) {
            minimumDistance = distances[0][i];
            farthest = i;
        }
    }
// Create a partial tour with two nodes
    tour[1] = farthest;
    visited_nodes[farthest] = true;
    currentSize++;
    tour[2] = 0;
// currentSize++;
    double *minimumAdditionalCosts = (double *) malloc(noOfThreads * sizeof(double));
    int *positions = (int *) malloc(noOfThreads * sizeof(int));
    int *nearestVertexes = (int *) malloc(noOfThreads * sizeof(int));
    int y = 0, threadID;
// Iterate through the rest of the nodes/coords
    while (currentSize < numOfCoords) {
        int max_index = 0;
        int position = 0;
        double max = DBL_MAX;
        double unvisitedmax = -1;
        int p;
        for (y = 0; y < noOfThreads; y++) {
            minimumAdditionalCosts[y] = 0;
            positions[y] = 0;
            nearestVertexes[y] = 0;
        }
/*
 * Step 3 - For each vertex Vi in the partial tour,
 * find an unvisited vertex Vp such that dist(Vi, Vp) is maximal
 */
#pragma omp parallel for collapse(2) private(i, p, threadID) shared(visited_nodes, distances, minimumAdditionalCosts, positions, nearestVertexes)
        for (i = 0; i < currentSize; i++) {
            for (p = 0; p < numOfCoords; p++) {
                threadID = omp_get_thread_num();
                if (!visited_nodes[p]) {
                    if (distances[tour[i]][p] > minimumAdditionalCosts[threadID]) {
                        minimumAdditionalCosts[threadID] = distances[tour[i]][p];
                        nearestVertexes[threadID] = p;
                    }
                }
            }
        }
        int x = 0;
        for (x = 0; x < noOfThreads; x++) {
            if (minimumAdditionalCosts[x] > unvisitedmax) {
                unvisitedmax = minimumAdditionalCosts[x];
                max_index = nearestVertexes[x];
            }
        }
/*
 * Step 4 : Insert Vp between two connected vertices in the partial tour Vi and Vi+1, where i is
 *  a position in the partial tour, such that dist(Vi, Vp) + dist(Vi+1, Vp) ?~H~R dist(Vi, Vi+1) is
 *  minimal.
 *
 */
        for (i = 0; i < currentSize; i++) {
            double additionalCost =
                    distances[tour[i]][max_index] + distances[tour[i + 1]][max_index] - distances[tour[i]][tour[i + 1]];
            if (additionalCost < max) {
                max = additionalCost;
                position = i + 1;
            }
        }
        currentSize++;
        visited_nodes[max_index] = true;
        for (i = currentSize; i > position; i--) {  //juhi
            tour[i] = tour[i - 1];
        }
        tour[position] = max_index;
        /*      for (i = 0; i < currentSize; i++) {
                  printf("visited array ########: %d \t %d\n", i, visited_nodes[i]);
              }*/
    }
// print the final tour
    for (i = 0; i <= currentSize; i++) {
        printf("%d \t", tour[i]);
    }
    int tourLength = i;
// Write the tour to the output file
    writeTourToFile(tour, tourLength, outputfile);
// Free allocated memory
/*    free(tour);
    free(visited_nodes);
free(minimumAdditionalCosts);
free(positions);
free(nearestVertexes);*/
}
