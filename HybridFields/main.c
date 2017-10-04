//
//  main.c
//  HybridFields
//
//  Created by Julian Manke on 21.01.16.
//  Copyright © 2016 Julian Manke. All rights reserved.
//


#include "Calculations.h"
#include "Particle.h"
#include "Grid.h"
#include "Fields.h"
#include <stdlib.h>
#include <stdio.h>
#include "Simulation.h"
#include <string.h>

int main() {
    
    ///calls testing method for the Liener-Wichert Fields. For this to work, change in "Calculations.c" within the method "writeElectricFieldToFile" the used value to "int zLW=70"h
    
    //Test if this method still works with the new xObserve[4]
    
    //LienerWiechert();
    
    ///calls testing method for the Maxwell-Pusher. For this to work, change in "Calculations.c" within the method "writeElectricFieldToFile" the used value to "int zMW" --> Chosse this value as input parameter?
    //maxwellPusher();
    
    //LienerWiechertPusher();
    hybridFieldsWithAnalyticPrecalculation();
   
    
    // Analytic Solution
    //LienerWiechertAnalyticSolution();
    
    
}
