#ifndef __ACC_HH__
#define __ACC_HH__

void updateGridNew(double grid[15][3], int cellNeighborMap[15][4]);
double* rk4StepWithNeighbors(double current[3], double neighbors[4][3], int neighborhoodSize);

#endif