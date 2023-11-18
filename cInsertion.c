#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>
#include<stdbool.h>


int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}

double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix) {

    // default private
    int i =0;
    int j = 0;
    for (i = 0; i < numOfCoords; i++) {
        for (j = 0; j < numOfCoords; j++) {

            double x1 = coordinates[i][0];
            double y1 = coordinates[i][1];
            double x2 = coordinates[j][0];
            double y2 = coordinates[j][1];

            double distance = calculateDistance(x1, y1, x2, y2);
            distanceMatrix[i][j] = distance;

        }

    }

    return distanceMatrix;
}

void cheapestInsertion(double **distanceMatrix, int numOfCoords, char *outputFileName)
{
    int visitedCount = 0;

    int *tour = (int*)malloc((numOfCoords+1)*sizeof(int));
    bool *visited = (bool*)malloc(numOfCoords*sizeof(bool));

    // Initialise with the first vertex
    tour[0] = 0;
    visited[0] = true;
    visitedCount++;

    // Find the nearest vertex
    double minimumDistance = DBL_MAX;

    int nearestVertex;
    int i = 0;
    for(i = 1 ; i <numOfCoords; i++)
    {
        if(distanceMatrix[0][i]< minimumDistance)
        {
            minimumDistance = distanceMatrix[0][i];
            nearestVertex = i;
        }
    }

    // Add the nearest vertex in the tour
    tour[1]= nearestVertex;
    visited[nearestVertex] = true;
    visitedCount++; // 2
    tour[2] = 0;

    while(visitedCount < numOfCoords)
    {
        double minimumAdditionalCost = DBL_MAX;

        int minN;
        int minUnvisited;
        // tour = {0,1,}
        for(i=0; i < visitedCount; i++)
        {
            // unvisited nodes
            int j = 0;
            for(j =0; j<numOfCoords; j++)
            {
                // check for unvisited nodes
                if(!visited[j])
                {
                    // j =2
                    double additionalCost = distanceMatrix[j][tour[i]]+ distanceMatrix[j][tour[i+1]] - distanceMatrix[tour[i]][tour[i+1]];
                    if(additionalCost < minimumAdditionalCost)
                    {
                        minimumAdditionalCost = additionalCost;
                        minN = i; // where to inset
                        minUnvisited = j; // what to insert
                    }
                }
            }
        }

        // Make space to add unvisited node to computed index
        for(i = visitedCount; i > minN; i--)
        {
            tour[i+1] = tour[i];
        }

        // add the node to tour
        tour[minN+1] = minUnvisited;
        visited[minUnvisited] = true;
        visitedCount++;

    }

    printf("Cheapest Insertion TSP Tour\n");

    double totalLength = numOfCoords+1;
    writeTourToFile(tour, totalLength, outputFileName);

}

int main(int argc, char *argv[]) {

    // Taking default file names if user didn't provide input
    char *fileName = "9_coords.coord";
    char *outputfileName = "output.txt";


    if (argc > 2) {
        fileName = argv[1];
        outputfileName = argv[2];
    }

    clock_t start, end;
    double time_taken;
    start = clock();


    int numOfCoords = readNumOfCoords(fileName);
    double **coordinates = readCoords(fileName, numOfCoords);

    printf("Number of coordinates:%d \n", numOfCoords);

    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));

    int i =0;
    for (i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);

    cheapestInsertion(distanceMatrix, numOfCoords, outputfileName);

    end = clock();
    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Cheapest insertion for %d nodes completed sequentially\n", numOfCoords);
    printf("The time taken is %fs .\n", time_taken);

    // Free memory
    for (i = 0; i < numOfCoords; i++) {
        free(coordinates[i]);
    }
    free(coordinates);

    for (i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);

    return 0;
}


