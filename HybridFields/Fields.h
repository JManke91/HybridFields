//
//  Fields.h
//  HybridFields
//
//  Created by Julian Manke on 12.01.17.
//  Copyright © 2017 Julian Manke. All rights reserved.
//

#ifndef Fields_h
#define Fields_h

#include <stdio.h>
#include "Grid.h"
#include "Simulation.h"

struct Fields {
    //double **electricField;
    //double **magneticField;
    //int numberofRowsElectricField;
    int long numberOfColumnsElectricField;
    int long numberOfColumnsAbsoluteElectricField;
    int long numberOfColumnsPLMCoefficients;
    double *electricField;//E-field
    double *magneticField;//H-field
    double *dField;
    double *bField;
    double *absoluteElectricField;
    //Parameters for UMPL (they get shifted differently for H and E fields according to "Tavfole"
    double *c1E;
    double *c2E;
    double *c3E;
    double *c4E;
    double *c5E;
    double *c6E;
    double *c1H;
    double *c2H;
    double *c3H;
    double *c4H;
    double *c5H;
    double *c6H;
    //Fields at box borders
    double **Hz_im1;
    double **Hy_im1;
    double **Hx_jm1;
    double **Hz_jm1;
    double **Hx_km1;
    double **Hy_km1;
    double **Ey_ip1;
    double **Ez_ip1;
    double **Ez_jp1;
    double **Ex_jp1;
    double **Ex_kp1;
    double **Ey_kp1;
};

typedef struct Fields Fields;
void initFields(Fields *Fields,Grid *Grid);
void allocateMemoryForFields(Fields *Fields);
void allocateMemoryForFieldsOnBoxBorders(Grid *Grid, Fields *Fields);
//void deallocareMemoryForFields(Fields *Fields);

#endif /* Fields_h */
