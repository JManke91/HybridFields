//
//  Grid.c
//  HybridFields
//
//  Created by Julian Manke on 05.01.17.
//  Copyright © 2017 Julian Manke. All rights reserved.
//

#include "Grid.h"
#include <stdlib.h>
#include <stdio.h>

void initGrid(Grid *Grid,double dx,double dy,double dz, int numberOfBoxesInX,int numberOfBoxesInY,int numberOfBoxesInZ,int numberOfGridPointsPerBoxInX,int numberOfGridPointsPerBoxInY,int numberOfGridPointsPerBoxInZ){
    allocateMemoryforXobserve(Grid);
    Grid->xobserve[0]=0; //x-component of x-observe
    Grid->xobserve[1]=0; //y-component of x-observe
    Grid->xobserve[2]=0; //z-component of x-observe
    Grid->dx=dx; //increment of Grid in x-direction
    Grid->dy=dy;
    Grid->dz=dz;
    Grid->numberOfGridPointsPerBoxInX = numberOfGridPointsPerBoxInX;
    Grid->numberOfGridPointsPerBoxInY = numberOfGridPointsPerBoxInY;
    Grid->numberOfGridPointsPerBoxInZ = numberOfGridPointsPerBoxInZ;
    Grid->numberOfBoxesInX = numberOfBoxesInX;
    Grid->numberOfBoxesInY = numberOfBoxesInY;
    Grid->numberOfBoxesInZ = numberOfBoxesInZ;
    Grid->Nx = numberOfBoxesInX*numberOfGridPointsPerBoxInX; //#of Grid points in x-direction
    Grid->Ny = numberOfBoxesInY*numberOfGridPointsPerBoxInY;
    Grid->Nz = numberOfBoxesInZ*numberOfGridPointsPerBoxInZ;
}

void allocateMemoryforXobserve(Grid *Grid) {
    Grid->xobserve=(double *)malloc(sizeof(double)*(4));
}

