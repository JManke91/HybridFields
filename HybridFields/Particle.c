//
//  Particle.c
//  HybridFields
//
//  Created by Julian Manke on 27.11.16.
//  Copyright © 2016 Julian Manke. All rights reserved.
//

#include "Particle.h"
#include "Calculations.h"
#include <stdlib.h>
#include <stdio.h>
//Complete Struct goes in as function variable and will be modified, with values which are plugged in
//explicit def. of function


void initParticle(particle *particle, double mass, int numberOfIterationSteps) {
    
    //particle->numberOfParticles = numberOfParticles;
    allocateMemory(particle,numberOfIterationSteps); //Allocate Memory for x,y,z Position Array of Particle
    //initValuesForArrays(particle, numberOfIterationSteps);
    particle->mass=mass;
    particle->iterationCount = 0;
    particle->historyLength = 0;
    particle->charge=1;
    particle->didParticleChangeBoxAfterPush = false;
    for(int i = 0; i<27; i++) {
        particle->boxIndicesOfNearFieldBoxesAfterPush[i] = -1;
        particle->boxIndicesOfNearFieldBoxesBeforePush[i] = -1;
    }
    for(int i = 0; i<4; i++) {
        particle->xRel[i] = 0;
        particle->uRel[i] = 0;
    }
    particle->oldBoxIndexBeforePush = -1;
    particle->newBoxIndexAfterPush = -1;
}


//// execute this on every array (yposition), etc... (which will be 2D now)
//void allocateMemoryfor2DArray(double*** array, int rows, int columns) {
//    *array = (double**)malloc(rows*sizeof(double*));
//    for(int i=0; i<rows; i++) {
//        (*array)[i] = (double*)malloc(columns*sizeof(double));
//    }
//}


///@brief allocate memory for history arrays for every single particle in simulation
void allocateMemory(particle *particle, int numberOfIterationSteps) {
    //    for(int i=1;i<=particle->numberOfParticles;i++) {
    particle->yposition = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->yposition == NULL) {
        printf("error allocating y array\n");
    }
    particle->xposition = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->xposition == NULL) {
        printf("error allocating x array\n");
    }
    particle->zposition = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->zposition == NULL) {
        printf("error allocating z array\n");
    }
    particle->xvelocity = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->xvelocity == NULL) {
        printf("error allocating x' array\n");
    }
    particle->yvelocity = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->yvelocity == NULL) {
        printf("error allocating y' array\n");
    }
    particle->zvelocity = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->zvelocity == NULL) {
        printf("error allocating z' array\n");
    }
    particle->velocityTime = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->velocityTime == NULL) {
        printf("error allocating velocityTime array\n");
    }
    particle->time = (double *) malloc(sizeof(double)*(numberOfIterationSteps + 1));
    if(particle->time == NULL) {
        printf("error allocating time array\n");
    }
    
    
//    particle->xRelHistory = (double **) malloc((numberOfIterationSteps + 1)* sizeof(double *));
//    if(particle->xRelHistory == NULL) {
//        printf("Error while trying to allocate memory for xRelHistory vector\n");
//    } else {
//        for(int i = 0; i < numberOfIterationSteps; i++) {
//            particle->xRelHistory[i] = (double *) malloc(4 * sizeof(double));
//            if(particle->xRelHistory[i] == NULL) {
//                printf("Error while trying to allocate memory for xRelHistory vector component\n");
//            }
//        }
//    }
    int i = 0;
    particle->xRelHistory = malloc((numberOfIterationSteps + 1)*sizeof(double *));
    if(particle->xRelHistory == NULL) {
        printf("error: out of memory.\n");
        return;
    }
    for(i = 0; i <= numberOfIterationSteps; i++) {
        particle->xRelHistory[i] = malloc((4)*sizeof(double));
        if(particle->xRelHistory[i] == NULL) {
            break;
        }
    }
    if(i != numberOfIterationSteps + 1) {
        printf("problem with allocation");
    }
    
    
    int j = 0;
    particle->uRelHistory = malloc((numberOfIterationSteps + 1)*sizeof(double *));
    if(particle->uRelHistory == NULL) {
        printf("error: out of memory.\n");
        return;
    }
    for(j = 0; j <= numberOfIterationSteps; j++) {
        particle->uRelHistory[j] = malloc((4)*sizeof(double));
        if(particle->uRelHistory[j] == NULL) {
            break;
        }
    }
    if(j != numberOfIterationSteps + 1) {
        printf("problem with allocation");
    }

    
//    particle->uRelHistory = (double **) malloc((numberOfIterationSteps + 1)* sizeof(double *));
//    if(particle->uRelHistory == NULL) {
//        printf("Error while trying to allocate memory for uRelHistory vector\n");
//    } else {
//        for(int i = 0; i < numberOfIterationSteps; i++) {
//            particle->uRelHistory[i] = (double *) malloc(4 * sizeof(double));
//            if(particle->uRelHistory[i] == NULL) {
//                printf("Error while trying to allocate memory for uRelHistory vector component\n");
//            }
//        }
//    }
    
}

void initValuesForArrays(particle *particle, int numberOfIterationSteps) {
    for(int i = 0; i < numberOfIterationSteps; i++) {
        for(int j = 0; j < 4; j++) {
            particle->xRelHistory[i][j] = 0;
            printf("xRelHistory[%i][%i] = %f\n", i, j, particle->xRelHistory[i][j]);
        }
    }
}


