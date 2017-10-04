//
//  Simulation.h
//  HybridFields
//
//  Created by Julian Manke on 13.01.17.
//  Copyright © 2017 Julian Manke. All rights reserved.
//

#ifndef Simulation_h
#define Simulation_h

#include <stdio.h>

//extern int number_of_Iterations_N_global;

//void simulation();
void LienerWiechert();
void maxwellPusher();
void LienerWiechertPusher();
void LienerWiechertAnalyticSolution();
void hybridFieldsWithAnalyticPrecalculation();
#endif /* Simulation_h */
