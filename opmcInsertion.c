#include <stdio.h>
#include "coordReader.c" // Include your custom header file
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>
#include<stdio.h>


double calculateDistance(double x1, double y1, double x2, double y2);
double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix);
void printDistanceMatrix(double **distanceMatrix, int numOfCoords);


double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}


void cheapestInsertion(double **distanceMatrix, int numOfCoords)
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
            // try testing omp parallel for
            minimumDistance = distanceMatrix[0][i];
            nearestVertex = i;
        }
    }

    // Add the nearest vertex in the tour
    tour[1]= nearestVertex;
    visited[nearestVertex] = true;
    visitedCount++; // 2
    tour[2] = 0;
    // paralell + test
    while(visitedCount < numOfCoords)
    {
        double minimumAdditionalCost = DBL_MAX;

        int minN;
        int minUnvisited;
        // tour = {0,1} paralel + collapse + crticial
        for(i=0; i < visitedCount; i++)
        {
            // unvisited nodes
            int j = 0;
            for(j =0; j<numOfCoords; j++)
            {
                // check for unvisited nodes
                if(!visited[j])
                {
                    // j =2 parallel + reduce
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

    printf("Cheapest Insertion TSP Tour By Cheap Tharki Harami Ketan:\n");


//    int supposedAnswer[] = {0, 11, 12, 3, 6, 10, 15, 1, 13, 5, 4, 2, 7, 8, 14, 9, 0};

    double totalLength = 0;
    double totalSupposedLength = 0;

    for ( i = 0; i <=numOfCoords; i++) {
        printf("%d ", tour[i]);
        if(i>0) {
            totalLength += distanceMatrix[tour[i]][tour[i - 1]];
        }
    }
    printf("\n");

    printf("%f", totalLength);
    printf("\n");

//    for ( i = 0; i <=numOfCoords; i++) {
//        printf("%d ", supposedAnswer[i]);
//        if(i>0) {
//            totalSupposedLength += distanceMatrix[supposedAnswer[i]][supposedAnswer[i - 1]];
//        }
//    }
    printf("\n");

//    printf("%f", totalSupposedLength);
    printf("\n");

}


int main(int argc, char *argv[]) {

    // taking default file name if user didn't provide input
    char *fileName = "9_coords.coord";

    if (argc > 1) {
        fileName = argv[1];
    }

    clock_t start, end;
    double time_taken;
    start = clock();

    printf("%s\n", fileName);


    int numOfCoords = readNumOfCoords("4096_coords.coord");
    double **coordinates = readCoords("4096_coords.coord", numOfCoords);

    printf("%dNunber of coords:", numOfCoords);
    printf("\n");

    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);

    cheapestInsertion(distanceMatrix, numOfCoords);

    end = clock();
    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("The time taken is %fs .\n", time_taken);

    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coordinates[i]);
    }
    free(coordinates);

    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);

    printf("Hello, World1223!\n");
    printf("%d", numOfCoords); // %d is the format specifier for integers
    return 0;
}



void printDistanceMatrix(double **distanceMatrix, int numOfCoords) {
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            printf("%f\t", distanceMatrix[i][j]);
        }
        printf("\n");
    }
}


double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix) {

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {

            double x1 = coordinates[i][0];
            double y1 = coordinates[i][1];
            double x2 = coordinates[j][0];
            double y2 = coordinates[j][1];


            double distance = calculateDistance(x1, y1, x2, y2);

            distanceMatrix[i][j] = distance;
            printf("%f\t", distanceMatrix[i][j]);

        }

        printf("\n");

    }

    return distanceMatrix;
}
