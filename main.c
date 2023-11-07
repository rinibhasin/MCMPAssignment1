#include <stdio.h>
#include "coordReader.h" // Include your custom header file
#include <float.h>
#include <stdlib.h>

double calculateDistance(double x1, double y1, double x2, double y2);
double **calculateDistanceMatrix(double **coordinates, int numOfCoords);
void printDistanceMatrix(double **distanceMatrix, int numOfCoords);
int findCheapestInsertion(double **distances, bool *visited, int *currentTour, int currentSize, int numOfCoords);
void cheapestInsertionTSP(double **distancematrix, int numOfCoords);



int main(int argc, char *argv[]) {

    char *filename = argv[1];
    printf("%s",filename);

    int numOfCoords = readNumOfCoords("9_coords.coord");
    double **coordinates = readCoords("9_coords.coord", numOfCoords);


    double **distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords);



    cheapestInsertionTSP(distanceMatrix,numOfCoords);
//    printDistanceMatrix(distanceMatrix, numOfCoords);
//    printf("\n");





    printf("Hello, World1223!\n");
    printf("%d", numOfCoords); // %d is the format specifier for integers
    return 0;
}






void printDistanceMatrix(double **distanceMatrix, int numOfCoords) {
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            printf("%.2f\t", distanceMatrix[i][j]);
        }
        printf("\n");
    }
}


double **calculateDistanceMatrix(double **coordinates, int numOfCoords) {
    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            double x1 = coordinates[i][0];
            double y1 = coordinates[i][1];
            double x2 = coordinates[j][0];
            double y2 = coordinates[j][1];

            double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

            distanceMatrix[i][j] = distance;
        }
    }

    return distanceMatrix;
}


// Function to find the index of the city to insert that minimizes tour length
int findCheapestInsertion(double **distances, bool *visited, int *currentTour, int currentSize, int numOfCoords) {
    int bestCity = -1;
    double bestIncrease = DBL_MAX;

    for (int city = 0; city < numOfCoords; city++) {
        if (!visited[city]) {
            for (int i = 0; i < currentSize; i++) {
                int next = (i + 1) % currentSize;
                double increase = distances[currentTour[i]][city] + distances[city][currentTour[next]] - distances[currentTour[i]][currentTour[next]];
                if (increase < bestIncrease) {
                    bestIncrease = increase;
                    bestCity = city;
                }
            }
        }
    }

    return bestCity;
}

// Function to solve the TSP using cheapest insertion
void cheapestInsertionTSP(double **distances, int numOfCoords) {

    int tour[numOfCoords];
    int currentSize = 1;
    bool visited[numOfCoords];

    // Start with the first city as the current tour
    tour[0] = 0;
    visited[0] = true;

    while (currentSize < numOfCoords) {
        int cityToInsert = findCheapestInsertion(distances, visited, tour, currentSize, numOfCoords);
        visited[cityToInsert] = true;

        // Find the position to insert the city in the current tour
        int insertPosition = -1;
        double bestIncrease = DBL_MAX;

        for (int i = 0; i < currentSize; i++) {
            int next = (i + 1) % currentSize;
            double increase = distances[tour[i]][cityToInsert] + distances[cityToInsert][tour[next]] - distances[tour[i]][tour[next]];

            if (increase < bestIncrease) {
                bestIncrease = increase;
                insertPosition = next;
            }
        }

        // Insert the city in the tour
        for (int i = currentSize; i > insertPosition; i--) {
            tour[i] = tour[i - 1];
        }

        tour[insertPosition] = cityToInsert;
        currentSize++;
    }

    // Print the final tour
    printf("Cheapest Insertion TSP Tour:\n");
    for (int i = 0; i < numOfCoords; i++) {
        printf("%d ", tour[i]);
    }
    printf("\n");
}

