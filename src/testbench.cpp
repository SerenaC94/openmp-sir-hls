#include "accelerator.hpp"
#include <stdio.h>

#ifdef __BAMBU__
#include <mdpi/mdpi_user.h>
#endif


int main() {
    double vectorS[15] = {0.999419,0.999517,0.913907,0.914528,0.998852,0.998493,0.895623,0.909938,0.894853,0.917034,0.93125,0.873875,0.916585,0.90144,0.904887};
    double vectorI[15] = {0,0,0.0416608,0.00325373,0,0,0.0110088,0.00185383,0.0889241,0.0313424,0.0640821,0.00110782,0.0146139,0.00717873,0.0918593};
    double vectorR[15] = {0.00058102,0.000482681,0.0444321,0.082218,0.00114778,0.00150748,0.0933686,0.0882079,0.0162229,0.0516231,0.00466816,0.125018,0.0688015,0.0913812,0.00325351};
    double grid[15][3];
    for (int i = 0; i < 15; i++)
    {
        grid[i][0] = vectorS[i];
        grid[i][1] = vectorI[i];
        grid[i][2] = vectorR[i];
    }

    int cellNeighborMap[15][4] = {
        {8,1,255,255},
        {9,0,2,255},
        {10,1,3,255},
        {11,2,4,255},
        {12,3,5,255},
        {13,4,6,255},
        {14,5,7,255},
        {15,6,8,255},
        {0,16,9,255},
        {1,17,8,10},
        {2,18,9,11},
        {3,19,10,12},
        {4,20,11,13},
        {5,21,12,14},
        {6,22,13,15}
    };

    printf("Grid:\n");
    for (int i = 0; i < 15; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            printf("%f, ", grid[i][j]);
        }
        printf("\n");
    }

    printf("Map:\n");
    for (int i = 0; i < 15; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            printf("%d, ", cellNeighborMap[i][j]);
        }
        printf("\n");
    }
    
    #ifdef __BAMBU_SIM__
    m_param_alloc(0, 15*3*sizeof(double));
    m_param_alloc(1, 15*4*sizeof(int));
    #endif
    
    // kernel to be accelerated
    updateGridNew(grid, cellNeighborMap);

    printf("Finished. New grid:\n");
    for (int i = 0; i < 15; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            printf("%f, ", grid[i][j]);
        }
        printf("\n");
    }

    return 0;
}