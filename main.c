#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>
#include<stdbool.h>

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

    printf("Farthest Insertion TSP Tour\n");


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

    // Can we make this allocation parallel
    int i =0;
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
