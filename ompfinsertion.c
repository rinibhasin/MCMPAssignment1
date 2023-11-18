#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>
#include<stdbool.h>
#include <omp.h>
#include <stdlib.h>


int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);


double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}

double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix) {

    int i, j;
    double x1, x2, y1, y2;
    double distance;

#pragma omp parallel for collapse(2) private(i, j, x1, y1, x2, y2, distance) shared(numOfCoords)
    for (i = 0; i < numOfCoords; i++) {
        for (j = 0; j < numOfCoords; j++) {

            x1 = coordinates[i][0];
            y1 = coordinates[i][1];
            x2 = coordinates[j][0];
            y2 = coordinates[j][1];

            distance = calculateDistance(x1, y1, x2, y2);
            distanceMatrix[i][j] = distance;

        }
    }

    return distanceMatrix;
}

void farthestInsertion(double **distanceMatrix, int numOfCoords, char *outputfilename)
{
    int visitedCount = 0;

    int *tour = (int*)malloc((numOfCoords+1)*sizeof(int));
    bool *visited = (bool*)malloc(numOfCoords*sizeof(bool));

    // Initialise with the first vertex
    tour[0] = 0;
    visited[0] = true;
    visitedCount++;

    // Find the farthest vertex
    double maximumDistance = DBL_MIN;

    int farthestVertex;
    int i = 0;
    for(i = 1 ; i <numOfCoords; i++)
    {
        if(distanceMatrix[0][i]> maximumDistance)
        {
            maximumDistance = distanceMatrix[0][i];
            farthestVertex = i;
        }
    }

    // Add the farthest vertex in the tour
    tour[1]= farthestVertex;
    visited[farthestVertex] = true;
    visitedCount++; // 2
    tour[2] = 0;

    int noOfThreads = omp_get_max_threads();

    double *farthestDistances = (double *) malloc(noOfThreads * sizeof(double));
    int *positions = (int *) malloc(noOfThreads * sizeof(int));
    int *farthestNodes = (int *) malloc(noOfThreads * sizeof(int));

    int y = 0, threadID;
    double unvisitedmax = -1;

    while(visitedCount < numOfCoords)
    {
        double farthestDistance = 0;
        int farthestNode;

        for (y = 0; y < noOfThreads; y++) {
            farthestDistances[y] = 0;
            positions[y] = 0;
            farthestNodes[y] = 0;
        }
        // tour = {0,1}
        int j = 0;
        #pragma omp parallel for collapse(2) private(i, j, threadID) shared(visited, distanceMatrix, farthestDistances, positions, farthestNodes)
        for(i=0; i <visitedCount; i++)
        {
            // Loop for nodes
            for(j=0; j<numOfCoords; j++)
            {
                threadID = omp_get_thread_num();
                // Checking for unvisited nodes
                if(!visited[j])
                {
                    // j =2
                    double currentDistance = distanceMatrix[j][tour[i]];
                    if(currentDistance > farthestDistances[threadID])
                    {
                        farthestDistances[threadID] = currentDistance;
                        farthestNodes[threadID] = j;
                    }
                }
            }
        }

        for (i = 0; i<noOfThreads; i++) {
            if (farthestDistances[i] > farthestDistance) {
                farthestDistance = farthestDistances[i];
                farthestNode = farthestNodes[i];
            }
        }

        double minimumAdditionalCost = DBL_MAX;
        int minN;
        for(i=0; i <visitedCount; i++)
        {
            // j =2
            double additionalCost = distanceMatrix[farthestNode][tour[i]]+ distanceMatrix[farthestNode][tour[i+1]] - distanceMatrix[tour[i]][tour[i+1]];
            if(additionalCost < minimumAdditionalCost)
            {
                minimumAdditionalCost = additionalCost;
                minN = i; // where to insert
            }
        }

        // Make space to add unvisited node to computed index
        for(i = visitedCount; i > minN; i--)
        {
            tour[i+1] = tour[i];
        }

        // Add the node to the tour
        tour[minN+1] = farthestNode;
        visited[farthestNode] = true;
        visitedCount++;

    }

    int tourLength = visitedCount+1;
    writeTourToFile(tour, tourLength, outputfilename);

}

int main(int argc, char *argv[]) {

    // taking default file names if user didn't provide input
    char *fileName = "16_coords.coord";
    char *outputfile = "output.txt";

    if (argc > 1) {
        fileName = argv[1];
        outputfile = argv[2];

    }

    double start, end;
    double time_taken;

    start = omp_get_wtime();;

    int numOfCoords = readNumOfCoords(fileName);
    double **coordinates = readCoords(fileName, numOfCoords);

    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));

    int i = 0;
    #pragma omp parallel for private(i)
    for (i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);
    farthestInsertion(distanceMatrix, numOfCoords, outputfile);

    end = omp_get_wtime();
    time_taken = (end - start);

    printf("Farthest insertion for %d nodes completed\n", numOfCoords);
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


