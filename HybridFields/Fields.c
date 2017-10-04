//
//  Fields.c
//  HybridFields
//
//  Created by Julian Manke on 12.01.17.
//  Copyright © 2017 Julian Manke. All rights reserved.
//

#include "Fields.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Grid.h"
#include "Calculations.h"

void initFields(Fields *Fields,Grid *Grid) {
    //Fields->numberofRowsElectricField=number_of_Iterations_N_global;
    Fields->numberOfColumnsElectricField=(Grid->Nx+1)*(Grid->Ny+1)*(Grid->Nz+1)*3;
    Fields->numberOfColumnsAbsoluteElectricField=(Grid->Nx+1)*(Grid->Ny+1)*(Grid->Nz+1);
    Fields->numberOfColumnsPLMCoefficients = Grid->Nx;
    allocateMemoryForFields(Fields);
    allocateMemoryForFieldsOnBoxBorders(Grid,Fields);
}

///@brief Memory Allocation for two dimensional electric field array. Structure of array (matrix):
///time =0=0th row {Ex,Ex,Ez(at pos. 000), Ex,Ey,Ez(at pos. 001),...,Ex,Ey,Ez(at pos.NxNyNz)} newline time=1=1strow {Ex,Ex,Ez(at pos. 000), Ex,Ey,Ez(at pos. 001),...,Ex,Ey,Ez(at pos.NxNyNz)} newline...etc. Equals: rows = N (#of time iterations) and columns = (Nx+1)*(Ny+1)*(Nz+1)*3(3 components for the field for each grid point including 0 --> +1)
void allocateMemoryForFields(Fields *Fields){
    
    
    //int rows=Fields->numberofRowsElectricField;
    int long columns=Fields->numberOfColumnsElectricField;
    int long absoluteColumns=Fields->numberOfColumnsAbsoluteElectricField;
    int long numberOfcolumnsPLMCoefficients = Fields->numberOfColumnsPLMCoefficients;
    
    Fields->electricField=(double *)malloc(sizeof(double)*(columns));
    Fields->magneticField=(double *)malloc(sizeof(double)*(columns));
    Fields->dField = (double *)malloc(sizeof(double)*(columns));
    Fields->bField = (double *)malloc(sizeof(double)*(columns));
    Fields->absoluteElectricField=(double *)malloc(sizeof(double)*(absoluteColumns));
    Fields->c1E = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c2E = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c3E = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c4E = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c5E = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c6E = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c1H = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c2H = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c3H = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c4H = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c5H = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    Fields->c6H = (double *)malloc(sizeof(double)*(numberOfcolumnsPLMCoefficients));
    
    //int long columnsAbsolute=Fields->numberOfColumnsAbsoluteElectricField;
    
    //malloc space for electric field
//    Fields->electricField=(double**)malloc(rows*sizeof(double*));
//    if(Fields->electricField==NULL) {printf("error assigning allocation for electricField rows\n");}
//    for(int i=0;i<rows;i++) {
//        Fields->electricField[i]=(double*)malloc(columns*sizeof(double*));
//        if(Fields->electricField[i]==NULL) {printf("error assigning electricField allocation in %i-th row for columns\n",i);}
//    }
    
    //malloc space for absoluteElectricFieldValues
//    Fields->absoluteElectricField=(double**)malloc(rows*sizeof(double*));
//    if(Fields->absoluteElectricField==NULL) {printf("error assigning allocation for absoluteElectric field for rows\n");}
//    for(int i=0;i<rows;i++) {
//        Fields->absoluteElectricField[i]=(double*)malloc(columnsAbsolute*sizeof(double*));
//        if(Fields->absoluteElectricField[i]==NULL) {printf("error assigning absoluteFied allocation in %i-th row for columns\n",i);}
//    }
}

///@brief allocates memory for the fields on the planes next to the box borders.
void allocateMemoryForFieldsOnBoxBorders(Grid *Grid, Fields *Fields) {
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    Fields->Hz_im1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Hy_im1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Hx_jm1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Hz_jm1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Hx_km1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Hy_km1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Ey_ip1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Ez_ip1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Ez_jp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Ex_jp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Ex_kp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    Fields->Ey_kp1 = (double**) malloc(numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ * sizeof(double*));
    
    if (Fields->Hz_im1 == NULL || Fields->Hy_im1 == NULL || Fields->Hx_jm1 == NULL || Fields->Hz_jm1 == NULL
        || Fields->Hx_km1 == NULL || Fields->Hy_km1 == NULL || Fields->Ey_ip1 == NULL || Fields->Ez_ip1 == NULL
        || Fields->Ez_jp1 == NULL || Fields->Ex_jp1 == NULL || Fields->Ex_kp1 == NULL || Fields->Ey_kp1 == NULL)
        printf("cannot allocate memory for the fields on box borders!\n");
    
    for (int i = 0; i < numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ; i++)
    {
        if ((Fields->Hz_im1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Hy_im1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Hx_jm1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Hz_jm1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Hx_km1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Hy_km1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Ey_ip1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Ez_ip1[i] = (double*) malloc(numberOfGridPointsForBoxInY * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Ez_jp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Ex_jp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInZ * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Ex_kp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        if ((Fields->Ey_kp1[i] = (double*) malloc(numberOfGridPointsForBoxInX * numberOfGridPointsForBoxInY * sizeof(double))) == NULL)
            printf("cannot allocate memory for the fields on box borders in for loop!\n");
        
    }
}

///@brief Memory Deallocation for two dimensional electric field array.
//void deallocareMemoryForFields(Fields *Fields) {
//    //deallocation of eletric field
//    for(int i=0;i<Fields->numberofRowsElectricField;i++){
//        free(Fields->electricField[i]);
//        Fields->electricField[i]=NULL;
//    }
//    free(Fields->electricField);
//    
//    //deallocation of absoluteElectricField
//    for(int i=0;i<Fields->numberofRowsElectricField;i++) {
//        free(Fields->absoluteElectricField[i]);
//        Fields->absoluteElectricField[i]=NULL;
//    }
//    free(Fields->absoluteElectricField);
    
    //deallocation of electric field array
    
//}
