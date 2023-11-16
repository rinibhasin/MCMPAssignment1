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



#define MIN(x, y) (((x) < (y)) ? (x) : (y))

double calculateDistance(double x1, double y1, double x2, double y2);
double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix);

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);


double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}


void cheapestInsertion(double **distanceMatrix, int numOfCoords)
{
    int visitedCount = 0;

    int *tour = (int*)malloc((numOfCoords+1)*sizeof(int));
    bool *visited = (bool*)malloc(numOfCoords*sizeof(bool));
    int m=0;
    for(m=0; m<numOfCoords; m++)
    {
        tour[m] = 0;
        visited[m] =0;
    }

    tour[numOfCoords] = 0;



    // Initialise with the first vertex
    tour[0] = 0;
    visited[0] = true;
    visitedCount++;

    // Find the nearest vertex
    double minimumDistance = DBL_MAX;

    int nearestVertex;
    int i = 0;

    int noOfThreads = 2;

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

    double *minimumAdditionalCosts = (double*)malloc(noOfThreads*sizeof(double));
    int *positions = (int*)malloc(noOfThreads*sizeof(int));
    int *nearestVertexes = (int*)malloc(noOfThreads*sizeof(int));

    int y =0;
    for(y =0;y< noOfThreads; y++)
    {
        minimumAdditionalCosts[y] = DBL_MAX;
        positions[y] =0;
        nearestVertexes[y]=0;
    }

    while(visitedCount < numOfCoords)
    {
        double minimumAdditionalCost = DBL_MAX;

        int minN;
        int minUnvisited;
        double additionalCost;
        int j;
        // tour = {0,1}

            #pragma omp parallel for collapse(2) private(i,j, additionalCost)
            for (i = 0; i < visitedCount; i++) {
                // unvisited nodes
                for (j = 0; j < numOfCoords; j++) {

                    int threadID = omp_get_thread_num();

                    printf("\n");
                    printf("%dthreadiD", threadID);
                    printf("\n");

                    // check for unvisited nodes
                    if (!visited[j]) {
                        // j =2
                        additionalCost = distanceMatrix[j][tour[i]] + distanceMatrix[j][tour[i + 1]] -
                                         distanceMatrix[tour[i]][tour[i + 1]];
                        if (additionalCost < minimumAdditionalCosts[threadID]) {
//                        minimumAdditionalCost = additionalCost;
//                        minN = i; // where to inset
//                        minUnvisited = j; // what to insert
                            printf("debug%f", additionalCost);
                            printf("\n");
                            minimumAdditionalCosts[threadID] = additionalCost;
                            printf("debug%f min additional cost stored", minimumAdditionalCosts[threadID]);
                            printf("\n");
                            positions[threadID] = i;
                            printf("debug%d******", positions[threadID]);
                            nearestVertexes[threadID] = j;
                        }
                    }
                }
            }



        int x=0;
        #pragma omp single
        for(x =0; x< noOfThreads; x++)
        {

            printf("The array is of length %d:", noOfThreads);
            printf("\n");
            printf("The element at %d:", x);
            printf("\n");
            printf("Minimum additional cost %f:", minimumAdditionalCosts[x]);
            printf("\n");
            printf("Minimum position %d:", positions[x]);
            printf("\n");
                printf("Minimum nearest vertex %d:", nearestVertexes[x]);

            if(minimumAdditionalCosts[x]< minimumAdditionalCost)
            {
                minimumAdditionalCost = minimumAdditionalCosts[x];
                minN = positions[x];
                minUnvisited = nearestVertexes[x];
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
    writeTourToFile(tour, totalLength, "output.txt");

}
void farthestInsertion(double **distanceMatrix, int numOfCoords)
{
    int visitedCount = 0;

    int *tour = (int*)malloc((numOfCoords+1)*sizeof(int));
    bool *visited = (bool*)malloc(numOfCoords*sizeof(bool));

    // Initialise with the first vertex
    tour[0] = 0;
    visited[0] = true;
    visitedCount++;

    // Find the nearest vertex
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

    // Add the nearest vertex in the tour
    tour[1]= farthestVertex;
    visited[farthestVertex] = true;
    visitedCount++; // 2
    tour[2] = 0;

    while(visitedCount < numOfCoords)
    {
        double farthestDistance = 0;
        int farthestNode;
        // tour = {0,1}
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
                    double currentDistance = distanceMatrix[j][tour[i]];
                    if(currentDistance > farthestDistance)
                    {
                        farthestDistance = currentDistance;
                        farthestNode = j; // where to inset
                    }
                }
            }
        }
        double minimumAdditionalCost = DBL_MAX;
        int minN;
        for(i=0; i < visitedCount; i++)
        {
            // j =2
            double additionalCost = distanceMatrix[farthestNode][tour[i]]+ distanceMatrix[farthestNode][tour[i+1]] - distanceMatrix[tour[i]][tour[i+1]];
            if(additionalCost < minimumAdditionalCost)
            {
                minimumAdditionalCost = additionalCost;
                minN = i; // where to inset
            }
        }

        // Make space to add unvisited node to computed index
        for(i = visitedCount; i > minN; i--)
        {
            tour[i+1] = tour[i];
        }
        printf("Current Visited Count %d\n", visitedCount);

        printf("Adding %d at position %d after %d:\n", farthestNode, minN + 1, tour[minN]);

        // add the node to tour
        tour[minN+1] = farthestNode;
        visited[farthestNode] = true;
        visitedCount++;

    }

    printf("Farthest Insertion TSP Tour \n");


    double totalLength = 0;

    for ( i = 0; i <=numOfCoords; i++) {
        printf("%d ", tour[i]);
        if(i>0) {
            totalLength += distanceMatrix[tour[i]][tour[i - 1]];
        }
    }
    printf("\n");

    printf("%f", totalLength);

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


    int numOfCoords = readNumOfCoords("9_coords.coord");
    double **coordinates = readCoords("9_coords.coord", numOfCoords);

    printf("%dNunber of coords:", numOfCoords);
    printf("\n");

    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));

    int i = 0;
    #pragma omp parallel for private(i)
    for (i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);

    cheapestInsertion(distanceMatrix, numOfCoords);
    //farthestInsertion(distanceMatrix, numOfCoords);

    end = clock();
    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
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

    printf("%d", numOfCoords); // %d is the format specifier for integers
    return 0;
}


double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix) {

   int i =0; int j =0;
   double x1 =0;
   double x2 =0;
   double y1 =0;
   double y2 =0;
   double distance = 0;

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
