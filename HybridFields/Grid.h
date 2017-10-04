//
//  Grid.h
//  HybridFields
//
//  Created by Julian Manke on 05.01.17.
//  Copyright © 2017 Julian Manke. All rights reserved.
//

#ifndef Grid_h
#define Grid_h

#include <stdio.h>

struct Grid {
    double dx;
    double dy;
    double dz;
    int Nx;
    int Ny;
    int Nz;
    int numberOfBoxesInX;
    int numberOfBoxesInY;
    int numberOfBoxesInZ;
    int numberOfGridPointsPerBoxInX;
    int numberOfGridPointsPerBoxInY;
    int numberOfGridPointsPerBoxInZ;
    //xobserve: struct array with 3 fixed entries, which will indicate the position at which retarded time (and therefore E and B fields) will be calculated
    double *xobserve;
};

typedef struct Grid Grid;
void initGrid(Grid *Grid,double dx,double dy,double dz, int numberOfBoxesInX,int numberOfBoxesInY,int numberOfBoxesInZ,int numberOfGridPointsPerBoxInX,int numberOfGridPointsPerBoxInY,int numberOfGridPointsPerBoxInZ);

void allocateMemoryforXobserve(Grid *Grid);



#endif /* Grid_h */
