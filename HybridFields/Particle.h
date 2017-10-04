//
//  Particle.h
//  HybridFields
//
//  Created by Julian Manke on 27.11.16.
//  Copyright © 2016 Julian Manke. All rights reserved.
//

#ifndef Particle_h
#define Particle_h
#include <stdio.h>
#include "Simulation.h"
#include <stdbool.h>

struct particle{
    double mass;
    //Needs to be 2d Array: First Index indicating particle number, second is access to array (like before)
    double *xposition;
    double *yposition;
    double *zposition;
    double *xvelocity;
    double *yvelocity;
    double *zvelocity;
    double *time;
    double xRel[4];
    double uRel[4];
    double **xRelHistory;
    double **uRelHistory;
    bool didParticleChangeBoxAfterPush;
    int boxIndicesOfNearFieldBoxesBeforePush[27];
    int boxIndicesOfNearFieldBoxesAfterPush[27];
    
    int oldBoxIndexBeforePush;
    int newBoxIndexAfterPush;
    // The first component of 4D (relativistic) velocity vector v^{\mu}=(\gamma,\gamma\vec{v})^{T} (c=0)
    double *velocityTime;
    //double *retardedTimeTest;
    double charge;
    // actual Phase of particle (i.e. real time values during simulation), including: actualPhase[0] = actual time, actualPhase[1] = current x position, actualPhase[2] = current y position, actualPhase[3] = current z position, actualPhase[4] = current x velocity, actualPhase[5] = current y velocity, actualPhase[6] = current z velocity
    double actualPhase[7];
    // New struct variable, indicating at how many entries above arrays have (equivalent to actual iteration step)- being updated in Nyström method
    int iterationCount;
    int historyLength;  
    // total # of particles in the simulation
    //int numberOfParticles;
    
};


typedef struct particle particle; //defines a new type "particle" and calls it "particle"
//Blueprint for the function defined in Particle.c

//struct particle particles[3]; //dynamic number of particles?

void initParticle(particle *particle, double mass, int numberOfIterationSteps);
void allocateMemory(particle *particle, int numberOfIterationSteps);
void allocateMemoryfor2DArray(double*** array, int rows, int columns);
void initValuesForArrays(particle *particle, int numberOfIterationSteps);

#endif /* Particle_h */
