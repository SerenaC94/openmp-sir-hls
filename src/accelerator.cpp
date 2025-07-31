#include "accelerator.hpp"

double* rk4StepWithNeighbors(double current[3], double neighbors[4][3], int neighborhoodSize) {
    double S = current[0];
    double I = current[1];
    double R = current[2];

    // Compute average neighbor infection level (Ij from neighbors)
    double totalI = 0.0;
    for (int i = 0; i < neighborhoodSize; i++) {
        totalI += neighbors[i][1];
    }
    double avgI = totalI / neighborhoodSize;

    // RK4 method (inlined versions of fS, fI, fR)
    // dt = 0.1, beta = 0.3, gammaRate = 0.1
    double k1_S = 0.1 * (-0.3 * S * avgI);
    double k1_I = 0.1 * (0.3 * S * avgI - 0.1 * I);
    double k1_R = 0.1 * (0.1 * I);

    double k2_S = 0.1 * (-0.3 * (S + 0.5 * k1_S) * avgI);
    double k2_I = 0.1 * (0.3 * (S + 0.5 * k1_S) * avgI - 0.1 * (I + 0.5 * k1_I));
    double k2_R = 0.1 * (0.1 * (I + 0.5 * k1_I));

    double k3_S = 0.1 * (-0.3 * (S + 0.5 * k2_S) * avgI);
    double k3_I = 0.1 * (0.3 * (S + 0.5 * k2_S) * avgI - 0.1 * (I + 0.5 * k2_I));
    double k3_R = 0.1 * (0.1 * (I + 0.5 * k2_I));

    double k4_S = 0.1 * (-0.3 * (S + k3_S) * avgI);
    double k4_I = 0.1 * (0.3 * (S + k3_S) * avgI - 0.1 * (I + k3_I));
    double k4_R = 0.1 * (0.1 * (I + k3_I));

    double newS = S + (k1_S + 2 * k2_S + 2 * k3_S + k4_S) / 6.0;
    double newI = I + (k1_I + 2 * k2_I + 2 * k3_I + k4_I) / 6.0;
    double newR = R + (k1_R + 2 * k2_R + 2 * k3_R + k4_R) / 6.0;

    static double new_cell[3];
    new_cell[0] = newS;
    new_cell[1] = newI;
    new_cell[2] = newR;

    return new_cell;
}

void updateGridNew(double grid[15][3], int cellNeighborMap[15][4]) {

    double tempGrid[15][3];
    #pragma omp parallel num_threads(THREAD_NUMBER)
    {
        #pragma omp for
        {
            for (int localIndex = 0; localIndex < 15; ++localIndex) {
                double currentCell[3];
                for(int i = 0; i < 3; i++)
                {
                    currentCell[i] = grid[localIndex][i];  
                } 
                double neighborsForUpdate[4][3];
                int neighborGlobalIds[4];
                for(int i = 0; i < 4; i++)
                {
                    neighborGlobalIds[i] = cellNeighborMap[localIndex][i];
                }
                int neighborhoodSize = 0;
                for (int neighborGlobalId = 0; neighborGlobalId < 4; ++neighborGlobalId) {
                    if (neighborGlobalIds[neighborGlobalId] < 16) {
                        neighborhoodSize++;
                        neighborsForUpdate[neighborGlobalId][0] = grid[neighborGlobalId][0];
                        neighborsForUpdate[neighborGlobalId][1] = grid[neighborGlobalId][1];
                        neighborsForUpdate[neighborGlobalId][2] = grid[neighborGlobalId][2];
                    }
                }
                double* temp;
                temp = rk4StepWithNeighbors(currentCell, neighborsForUpdate, neighborhoodSize);
                tempGrid[localIndex][0] = temp[0];
                tempGrid[localIndex][1] = temp[1];
                tempGrid[localIndex][2] = temp[2];
            }
        }
        #pragma omp for
            for (int localIndex = 0; localIndex < 15; ++localIndex) {
                grid[localIndex][0] = tempGrid[localIndex][0];
                grid[localIndex][1] = tempGrid[localIndex][1];
                grid[localIndex][2] = tempGrid[localIndex][2];
        }
    }
        
    // Apply normalization as part of the computation step
    for (int i = 0; i < 15; i++) {
        double sum_sir = grid[i][0] + grid[i][1] + grid[i][2];
        if (sum_sir > 1e-9) {
            grid[i][0] = (grid[i][0] / sum_sir);
            grid[i][1] = (grid[i][1] / sum_sir);
            grid[i][2] = (grid[i][2] / sum_sir);
        } else {
            grid[i][0] = (1.0);
            grid[i][1] = (0.0);
            grid[i][2] = (0.0);
        }
    }

}