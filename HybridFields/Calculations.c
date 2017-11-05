//
//  Calculations.c
//  HybridFields
//
//  Created by Julian Manke on 27.11.16.
//  Copyright © 2016 Julian Manke. All rights reserved.
//

#include "Calculations.h"
#include <stdbool.h>
#include "Particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Grid.h"
#include "Fields.h"
#include <string.h>

///@brief maxwellPusher for H field.
///@remark Method only used for test purposes. H field is calculated via positive curl. Therefore value of E on the right side of H on the grid is required. Thus stop i,j,k with at n - 1 where n denotes the numberOfGridPoints.
///@param Grid pointer to Grid struct
///@param dt time increment
void pushHFieldOnGrid(Fields *Fields, Grid *Grid,double dt) {
    double Ex_ijk;
    double Ey_ijk;
    double Ez_ijk;
    double Ey_ijkp1;
    double Ez_ijp1k;
    double Ez_ip1jk;
    double Ex_ijkp1;
    double Ex_ijp1k;
    double Ey_ip1jk;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    
    int nx = Grid->Nx;
    int ny = Grid->Ny;
    int nz = Grid->Nz;
    
    int i, j, k;
    for (i = 1; i < nx - 1; i++) {
        for (j = 1; j < ny - 1; j++) {
            for (k = 1; k < nz - 1; k++) {
                Ey_ijkp1 = Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k + 1) + 1];
                Ey_ijk = Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 1];
                Ez_ijp1k = Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j + 1) + 3 * (k) + 2];
                Ez_ijk = Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 2];
                Ez_ip1jk = Fields->electricField[3 * (nz) * (ny) * (i + 1) + 3 * (nz) * (j) + 3 * (k) + 2];
                Ex_ijkp1 = Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k + 1) + 0];
                Ex_ijk = Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 0];
                Ex_ijp1k = Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j + 1) + 3 * (k) + 0];
                Ey_ip1jk = Fields->electricField[3 * (nz) * (ny) * (i + 1) + 3 * (nz) * (j) + 3 * (k) + 1];
                
                
                double add_x = cnz * (Ey_ijkp1 - Ey_ijk) - cny * (Ez_ijp1k - Ez_ijk);
                double add_y = cnx * (Ez_ip1jk - Ez_ijk) - cnz * (Ex_ijkp1 - Ex_ijk);
                double add_z = cny * (Ex_ijp1k - Ex_ijk) - cnx * (Ey_ip1jk - Ey_ijk);
                
                Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 0] += add_x;
                Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 1] += add_y;
                Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 2] += add_z;
            }
        }
    }
}


///@brief maxwellPusher for E field. This method is only used for testing purposes. See pushEFieldInsideBoxes for the final version. Fields will be calculated going through the grid xyz. Actual implementation may require to go through grid like zxy
///@remark Method only used for test purposes. E field is calculated via negative curl. Therefore value of B on the left side of E on the grid is required. Thus start i,j,k with 1
///@param Grid pointer to Grid struct
///@param dt time increment to calculate electric field for
void pushEFieldOnGrid(Fields *Fields,Grid *Grid, double dt){
    double Hx_ijk;
    double Hy_ijk;
    double Hz_ijk;
    double Hz_ijm1k;
    double Hy_ijkm1;
    double Hx_ijkm1;
    double Hz_im1jk;
    double Hy_im1jk;
    double Hx_ijm1k;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    
    int nx = Grid->Nx;
    int ny = Grid->Ny;
    int nz = Grid->Nz;
    
    int i, j, k;
    for (i = 1; i < nx - 1; i++) {
        for (j = 1; j < ny - 1; j++) {
            for (k = 1; k < nz - 1; k++) {
                Hx_ijk = Fields->magneticField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * k + 0];
                Hy_ijk = Fields->magneticField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * k + 1];
                Hz_ijk = Fields->magneticField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * k + 2];
                
                Hz_ijm1k = Fields->magneticField[3 * (nz) * (ny) * i + 3 * (nz) * (j - 1) + 3 * k + 2];
                Hy_ijkm1 = Fields->magneticField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (k - 1) + 1];
                Hx_ijkm1 = Fields->magneticField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (k - 1) + 0];
                Hz_im1jk = Fields->magneticField[3 * (nz) * (ny) * (i - 1) + 3 * (nz) * j + 3 * k + 2];
                Hy_im1jk = Fields->magneticField[3 * (nz) * (ny) * (i - 1) + 3 * (nz) * (j) + 3 * (k) + 1];
                Hx_ijm1k = Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j - 1) + 3 * (k) + 0];
                
                
                double val_x = cny * (Hz_ijk - Hz_ijm1k) - cnz * (Hy_ijk - Hy_ijkm1);
                double val_y = cnz * (Hx_ijk - Hx_ijkm1) - cnx * (Hz_ijk - Hz_im1jk);
                double val_z = cnx * (Hy_ijk - Hy_im1jk) - cny * (Hx_ijk - Hx_ijm1k);
                
                Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 0] += val_x;
                Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 1] += val_y;
                Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 2] += val_z;
            }
        }
    }
}


///@brief Checks, if a value is contained in an Array
///@param value Value to be checked
///@param arr Array, in which function is looking for value
///@size Size of the Array
bool doesArrayContainValue(int value, int arr[], int size) {
    for(int i=0;i<size;i++) {
        if(arr[i] == value) {
            return true;
        }
    }
    return false;
}


///@brief Sets all LW fields within a whole box, given by index "boxIndex" to zero.
void setLWFieldsInBoxToZero(Grid *Grid, Fields *Fields, int boxIndex) {
    double xObserve[4];
    // Def. first component of xObserve for convenience (actual value is not needed for calculation)
    xObserve[0]=0;
    // first, from box index "boxindex", calculate ib,jb,kb: int ib(index) = boxIndex/(#BoxesZ*#BoxesY), int jb(boxIndex) = boxIndex/#BoxesY % BoxesY, int kb = boxIndex % #BoxesZ
    int ib = boxIndex/(Grid->numberOfBoxesInZ*Grid->numberOfBoxesInY);
    int jb = boxIndex/(Grid->numberOfBoxesInY) % (Grid->numberOfBoxesInY);
    int kb = boxIndex % (Grid->numberOfBoxesInZ);
    // from ib,jb,kb, calculate corresponding indices on grid at the bottom left corner of box: int ic = ib*#GridPointsPerBoxInX,...etc.
    int i = ib*Grid->numberOfGridPointsPerBoxInX;
    int j = jb*Grid->numberOfGridPointsPerBoxInY;
    int k = kb*Grid->numberOfGridPointsPerBoxInZ;
    // Define Parameters until which Loops should go
    int imax = i+Grid->numberOfGridPointsPerBoxInX;
    int jmax = j+Grid->numberOfGridPointsPerBoxInY;
    int kmax = k+Grid->numberOfGridPointsPerBoxInZ;
    // Starting at ic, jc, kc loop through a whole box adding fields at each grid point
    // ATTENTION: Check, if indices ic, jc, kc are set correctly
    for(int ic = i;ic<imax;ic++) {
        for(int jc = j;jc<jmax;jc++) {
            for(int kc = k;kc<kmax;kc++) {
                xObserve[1] = ic*Grid->dx;
                xObserve[2] = jc*Grid->dy;
                xObserve[3] = kc*Grid->dz;
                
                Fields->electricField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 0] = 0;
                Fields->electricField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 1] = 0;
                Fields->electricField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 2] = 0;
                
                Fields->magneticField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 0] = 0;
                Fields->magneticField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 1] = 0;
                Fields->magneticField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 2] = 0;
            }
        }
    }
}


///@brief Subs Linard Wiechert Fields within a whole box, given by index "boxIndex"
void subLWFieldsInBox(Grid *Grid, Fields *Fields, particle *particle, int boxIndex, double time) {
    // Def. xObserve (that's, where fields are calculated)
    //printf("subtrating LW in box %d\n", boxIndex);
    double xObserve[4];
    // Def. first component of xObserve for convenience
    xObserve[0]=time;
    // first, from box index "boxindex", calculate ib,jb,kb: int ib(index) = boxIndex/(#BoxesZ*#BoxesY), int jb(boxIndex) = boxIndex/#BoxesY % BoxesY, int kb = boxIndex % #BoxesZ
    int ib = boxIndex/(Grid->numberOfBoxesInZ*Grid->numberOfBoxesInY);
    int jb = boxIndex/(Grid->numberOfBoxesInY) % (Grid->numberOfBoxesInY);
    int kb = boxIndex % (Grid->numberOfBoxesInZ);
    // from ib,jb,kb, calculate corresponding indixes on grid at the bottom left corner of box: int ic = ib*#GridPointsPerBoxInX,...etc.
    int i = ib*Grid->numberOfGridPointsPerBoxInX;
    int j = jb*Grid->numberOfGridPointsPerBoxInY;
    int k = kb*Grid->numberOfGridPointsPerBoxInZ;
    // Define Parameters until which Loops should go
    int imax = i+Grid->numberOfGridPointsPerBoxInX;
    int jmax = j+Grid->numberOfGridPointsPerBoxInY;
    int kmax = k+Grid->numberOfGridPointsPerBoxInZ;
    // Starting at ic, jc, kc loop through a whole box adding fields at each grid point
    // ATTENTION: Chef, if indices ic, jc, kc are set correctly
    for(int ic = i;ic<imax;ic++) {
        for(int jc = j;jc<jmax;jc++) {
            for(int kc = k;kc<kmax;kc++) {
                xObserve[1] = ic*Grid->dx;
                xObserve[2] = jc*Grid->dy;
                xObserve[3] = kc*Grid->dz;
                subLWField(Grid, particle, &Fields->electricField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 0], xObserve, 0);
                subLWField(Grid, particle, &Fields->electricField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 1], xObserve, 1);
                subLWField(Grid, particle, &Fields->electricField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 2], xObserve, 2);
                subLWField(Grid, particle, &Fields->magneticField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 0], xObserve, 3);
                subLWField(Grid, particle, &Fields->magneticField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 1], xObserve, 4);
                subLWField(Grid, particle, &Fields->magneticField[3 * (Grid->Nz) * (Grid->Ny) * (ic) + 3 * (Grid->Nz) * (jc) + 3 * kc + 2], xObserve, 5);
            }
        }
    }
}


///@brief Constructs a sample pluse for Maxwell Pusher (multidimensional Gaussian pulse) onto the grid for testing purposes.
///@param Grid pointer to Grid struct
void initSamplePulseOnGridGausssian(Grid *Grid, Fields *Fields) {
    
    int nx = Grid->Nx;
    int ny = Grid->Ny;
    int nz = Grid->Nz;
    double mu = 64; // peak of Gaussian pulse
    double sigma = 10; //wideness of Gaussian pulse --> 5 is a good choice
    double w = 15.7; // Frequency of modulated Sine pulse. Border k vector (at which group velocity is supposed to be zero): k_b = pi/\Delta x \approx 15,7
    
    for (int i = 0; i < nx; i++){ // 32*2 and 32*6 // 20 to 108
        for (int j = 0; j < ny; j++){
            for (int k = 0; k < nz; k++){
                Fields->electricField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2] =  exp(-(pow((i - mu) * Grid->dx, 2)/sigma + pow((j - mu) * Grid->dy, 2)/sigma + pow((k - mu) * Grid->dz, 2)/sigma)) * sin(w * i);
                Fields->magneticField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1] =  exp(-(pow((i - mu) * Grid->dx, 2)/sigma + pow((j - mu) * Grid->dy, 2)/sigma + pow((k - mu) * Grid->dz, 2)/sigma)) * sin(w * i);
            }
        }
    }
}


void initializeParticle(particle *particle, double t0, double x0, double y0, double z0, double xDot0, double yDot0, double zDot0) {
    
    particle->xRel[0] = t0;
    particle->xRel[1] = x0;
    particle->xRel[2] = y0;
    particle->xRel[3] = z0;
    
    particle->uRel[0] = sqrt(1+pow(xDot0,2)+pow(yDot0,2)+pow(zDot0,2));
    particle->uRel[1] = xDot0;
    particle->uRel[2] = yDot0;
    particle->uRel[3] = zDot0;
    
    particle->xRelHistory[0][0] = t0;
    particle->xRelHistory[0][1] = x0;
    particle->xRelHistory[0][2] = y0;
    particle->xRelHistory[0][3] = z0;
    
    particle->uRelHistory[0][0] = particle->uRel[0];
    particle->uRelHistory[0][1] = xDot0;
    particle->uRelHistory[0][2] = yDot0;
    particle->uRelHistory[0][3] = zDot0;
    
    particle->iterationCount = 0;
    particle->historyLength = 0;
}

bool double_equals(double a, double b) {
    double epsilon = 0.001;
    return fabs(a-b) < epsilon;
}


///@brief takes electric field (1D) vector and writes it into text file, formatted for constant z value in x,y plane, i.e. {(00,10,20,...),(01,11,21,...),...}, for (...) being one matrix row. Writing: z,x,y
///@param particle Used for the constant (absolute) z-coordinate, where particle is supposed to be plotted, e.g. if particle starts at (1,1,1.5) and circles around in xy plane, the zCoordinate is just 1.5
void writeElectricFieldToFile(Grid *Grid, particle *particle, Fields *Fields, int index) {
    //Create correct Filename for each time iteration
    char filename[30]="electricFieldAtTime";
    char character[20];
    //int nx = Grid->Nx;
    int ny = Grid->Ny;
    int nz = Grid->Nz;
    
    //Convert input String "index" to character
    sprintf(character,"%d",index);
    //Combine "filename" and "index"
    strcat(filename,character);
    //Combine new filename and ".txt"
    strcat(filename,".txt"); // Now "filename" = "electricFieldAtTimei.txt", where i is actual time index
    int zLW = (particle->xRel[3])/(Grid->dz);
     //int zComponentPuse = 64; // for pulse propagation (Maxwell-Pusher method)
    
    FILE *file=fopen(filename,"w");
    
    // next row
    for(int j=0; j<Grid->Ny; j++){
        //fills columns for one row
        for(int i=0;i<Grid->Nx;i++) {
            double Ex = Fields->electricField[3 * nz * ny * i + 3 * nz * j + 3 * zLW + 0];
            double Ey = Fields->electricField[3 * nz * ny * i + 3 * nz * j + 3 * zLW + 1];
            double Ez = Fields->electricField[3 * nz * ny * i + 3 * nz * j + 3 * zLW + 2];
            double absoluteValueOfField = pow(Ex,2)+pow(Ey,2)+pow(Ez,2);
            fprintf(file,"%f\t",absoluteValueOfField);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}


///@brief loops through the entire E and B array and writes Bx,By,Bz and (or) Ex,Ey,Ez to seperate files. For Fourier Analysis, not the intensity \propto  |E^2|, etc 
///@param filename pointer to a char. Gets modified inside the method
///@param index outer loop index. Is used to name the output file
///@param plotE set to true, if you want to write E-field to file
///@param plotB set to true, if you want to write B-field to file
///@throws ERROR: Could not open file for E or B field
void writeFieldComponentsForFourierAnalysisToFile(Grid *Grid, Fields *Fields, char *filename, int index, int planeForPlotting, bool plotE, bool plotB) {
    printf("Writing field components for Fourier Analysis to file ...\n");
    FILE *fid_Ex = NULL;
    FILE *fid_Ey = NULL;
    FILE *fid_Ez = NULL;
    FILE *fid_Hx = NULL;
    FILE *fid_Hy = NULL;
    FILE *fid_Hz = NULL;
    
    if (plotE){
        sprintf(filename, "E_field_x%d", index);
        strcat(filename, ".txt");
        fid_Ex = fopen(filename,"w");
        if (fid_Ex == NULL){
            printf("ERROR: Could not open file E_field_x!");
        }
        sprintf(filename, "E_field_y%d", index);
        strcat(filename, ".txt");
        fid_Ey = fopen(filename,"w");
        if (fid_Ey == NULL){
            printf("ERROR: Could not open file E_field_y!");
        }
        sprintf(filename, "E_field_z%d", index);
        strcat(filename, ".txt");
        fid_Ez = fopen(filename,"w");
        if (fid_Ez == NULL){
            printf("ERROR: Could not open file E_field_z!");
        }
    }
    if(plotB){
        sprintf(filename, "B_field_x%d", index);
        strcat(filename, ".txt");
        fid_Hx = fopen(filename,"w");
        if (fid_Hx == NULL){
            printf("ERROR: Could not open file B_field_x!");
        }
        sprintf(filename, "B_field_y%d", index);
        strcat(filename, ".txt");
        fid_Hy = fopen(filename,"w");
        if (fid_Hy == NULL){
            printf("ERROR: Could not open file B_field_y!");
        }
        sprintf(filename, "B_field_z%d", index);
        strcat(filename, ".txt");
        fid_Hz = fopen(filename,"w");
        if (fid_Hz == NULL){
            printf("ERROR: Could not open file B_field_z!");
        }
    }
    
    else{
        int nx = Grid->Nx;
        int ny = Grid->Ny;
        int nz = Grid->Nz;
        int k = planeForPlotting;
        double Ex, Ey, Ez, Hx, Hy, Hz;
        
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                Hx = Fields->magneticField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0];
                Hy = Fields->magneticField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1];
                Hz = Fields->magneticField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2];
                
                Ex = Fields->electricField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 0];
                Ey = Fields->electricField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 1];
                Ez = Fields->electricField[3 * nz * ny * (i) + 3 * nz * (j) + 3 * (k) + 2];
                if(plotE){
                    fprintf(fid_Ex, "%.18f\t", Ex);
                    fprintf(fid_Ey, "%.18f\t", Ey);
                    fprintf(fid_Ez, "%.18f\t", Ez);
                }
                if (plotB){
                    fprintf(fid_Hx, "%.18f\t", Hx);
                    fprintf(fid_Hy, "%.18f\t", Hy);
                    fprintf(fid_Hz, "%.18f\t", Hz);
                }
            }
            if(plotE){
                fprintf(fid_Ex,"\n");
                fprintf(fid_Ey,"\n");
                fprintf(fid_Ez,"\n");
            }
            if(plotB){
                fprintf(fid_Hx,"\n");
                fprintf(fid_Hy,"\n");
                fprintf(fid_Hz,"\n");
            }
        }
    }
    fclose(fid_Ex);
    fclose(fid_Ey);
    fclose(fid_Ez);
    fclose(fid_Hx);
    fclose(fid_Hy);
    fclose(fid_Hz);
}

///@brief writing detailed simulation info (for grid and all particle positions & velocities) into file, such that simulation is easily to be identified later on.
void writeDetailedSimulationInfoToFile(particle *particles, Grid *Grid, int numberOfParticles, double precalculationTime, double calculationTime) {
    FILE *detailedInfo = fopen("detailedSimulationInfo.txt", "w");
    // general info
    fprintf(detailedInfo, "Number of particles = %i\nPrecalculation time = %f\nCalculation time = %f\n\ndx = %f\ndy = %f\ndz = %f\nNumber of grid points per box in x = %d\nNumber of grid points per box in y = %d\nNumber of grid points per box in z = %d\nNumber of boxes in x = %d\nNumber of boxes in y = %d\nNumber of boxes in z = %d\n\n", numberOfParticles, precalculationTime, calculationTime, Grid->dx, Grid->dy, Grid->dz, Grid->numberOfGridPointsPerBoxInX, Grid->numberOfGridPointsPerBoxInY, Grid->numberOfGridPointsPerBoxInZ, Grid->numberOfBoxesInX, Grid->numberOfBoxesInY, Grid->numberOfBoxesInZ);
    // particle info
    for(int p = 0; p < numberOfParticles; p++) {
        fprintf(detailedInfo, "Particle %i initial positions...x[1] = %f, x[2] = %f, x[3] = %f\nParticle %i initial velocities...u[1] = %f, u[2] = %f, u[3] = %f\n", p, particles[p].xRel[1], particles[p].xRel[2], particles[3].xRel[3], p, particles[p].uRel[1], particles[p].uRel[2], particles[p].uRel[3]);
    }
    fclose(detailedInfo);
}


///@brief calculates the Minkowski Product of two 4D Minkowski Vectors with sign convention (+,-,-,-)
double minkowskiProduct(double x[4], double y[4]) {
    double product = x[0]*y[0]-x[1]*y[1]-x[2]*y[2]-x[3]*y[3];
    return product;
}


///@brief calculated the standard scalar product for a 4D vector
double standard4DScalarProduct(double inputOne[4], double inputTwo[4]) {
    double product = inputOne[0] * inputTwo[0] + inputOne[1] * inputTwo[1] + inputOne[2] * inputTwo[2] + inputOne[3] * inputTwo[3];
    return product;
}

void freeMemoryforArray(double array[]) {
    free(array);
}

///@brief calculates all neccessary parameters to be able to calculate the Liener Wiechert Fields
void calculateLienerWiechertParameters(double intersectionPoint[4], double velocityAtIntersectionPoint[4],double xObserve[4],double n[3],double *R,double beta[3],double *gamma_sq) {
    //Calculation of \vec{n}_{pq}=(\vec{x}_{observe}(t_{observe})-\vec{x}_{particle}(t_{retarded}))
    n[0]=xObserve[1]-intersectionPoint[1];
    n[1]=xObserve[2]-intersectionPoint[2];
    n[2]=xObserve[3]-intersectionPoint[3];
    
    *R=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    
    n[0]/=*R;
    n[1]/=*R;
    n[2]/=*R;
    
    *gamma_sq=pow(velocityAtIntersectionPoint[0],2);
    
    //beta
    beta[0]=(velocityAtIntersectionPoint[1])/(velocityAtIntersectionPoint[0]);
    beta[1]=(velocityAtIntersectionPoint[2])/(velocityAtIntersectionPoint[0]);
    beta[2]=(velocityAtIntersectionPoint[3])/(velocityAtIntersectionPoint[0]);
}

///@brief calculates all neccessary parameters to be able to calculate the Liener Wiechert Fields
void calculateLienerWiechertParametersForCertainComponent(double intersectionPoint[4], double velocityAtIntersectionPoint[4],double xObserve[4],double n[3],double *R,double beta[3],double *gamma_sq) {
    //Calculation of \vec{n}_{pq}=(\vec{x}_{observe}(t_{observe})-\vec{x}_{particle}(t_{retarded}))
    n[0]=xObserve[1]-intersectionPoint[1];
    n[1]=xObserve[2]-intersectionPoint[2];
    n[2]=xObserve[3]-intersectionPoint[3];
    
    *R=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    
    n[0]/=*R;
    n[1]/=*R;
    n[2]/=*R;
    
    *gamma_sq=pow(velocityAtIntersectionPoint[0],2);
    
    //beta
    beta[0]=(velocityAtIntersectionPoint[1])/(velocityAtIntersectionPoint[0]);
    beta[1]=(velocityAtIntersectionPoint[2])/(velocityAtIntersectionPoint[0]);
    beta[2]=(velocityAtIntersectionPoint[3])/(velocityAtIntersectionPoint[0]);
}

///@brief calculates the cross product of two vectors "a" and "b" and saves the result in array "result"
///@param a first vector of cross product
///@param b second vector of cross product
///@param result the result vector which  is automatically available in outer scope
void crossProduct(double a[3], double b[3], double result[3]){
    result[0] = a[1]*b[2]-a[2]*b[1];
    result[1] = a[2]*b[0]-a[0]*b[2];
    result[2] = a[0]*b[1]-a[1]*b[0];
}

///@brief In the context of analytical precalulation of fields, after reversed particle push, the particle push in "desired" direction starts. But originally not at desired initial conditions, but at positions obtained from reversed particle push. Fronm here, particle gets pushed via Nystrom up to the desired predefined initial conditions. The particle history up to those predefined initial positions is taken to calculate LW fields analytically on the whole grid (only in the very last time step, because right after that, the hybrid field calculation with Maxwell Pushers set in!).
void precalculateFieldsForGivenPrecalculationTime(particle *particles, Grid *Grid, Fields *Fields, int numberOfParticles, int numberOfPrecalculationSteps, FILE *rect, double dt, double tN, double precalculationTime, double Bext[3], double Eext[3]) {
    
    // predefine stuff
//    int numberOfBoxesInX = Grid->numberOfBoxesInX;
//    int numberOfBoxesInY = Grid->numberOfBoxesInY;
//    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;

    for(int step = 0; step < precalculationTime / dt; step++) {
        //*timeIterationStep += 1;
        
        printf("Precalculation of time step %i of %i\n", step, numberOfPrecalculationSteps);
        
        for(int h=0; h < numberOfParticles; h++) {
            
            // write casted particle position (for rect) into file ============= --> ADJUST FOR MORE THAN TWO PARTICLES
            int numberBoxesX = (int) (particles[h].xRel[1]/(Grid->numberOfGridPointsPerBoxInX * Grid->dx));
            int numberBoxesY = (int) (particles[h].xRel[2]/(Grid->numberOfGridPointsPerBoxInY * Grid->dy));
            
            double cornerX = numberBoxesX * (Grid->numberOfGridPointsPerBoxInX * Grid->dx);
            double cornerY = numberBoxesY * (Grid->numberOfGridPointsPerBoxInY * Grid->dy);
            fprintf(rect, "%f %f\n", cornerX, cornerY);
            // =================================================================
            addCurrentStateToParticleHistory(&particles[h], step);
            
            //double tau = dt/particles[h].uRel[0];
            //Nystrom(particles, Fields, Grid, step, tau, numberOfParticles, h, dt, tN, Bext, Eext, true);
            
            updateVelocityWithBorisPusher(particles, Grid, Fields, numberOfParticles, h, Eext, Bext, dt, step);
            updateLocation(particles, Grid, dt, h, step);
            
            
            // just at the very last time step, calculate fields for whole particle trajectory
//            if(double_equals(step * dt, (precalculationTime - dt))) {
//                // code
//                for(int boxIndex = 0; boxIndex < numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ; boxIndex++) {
//                    if(!boxIsInNearFieldOfParticle(Grid, &particles[h], boxIndex)) {
//                        double evaluationTime = step * dt;
//                        addLWFieldsInBoxforTest(Grid, Fields, &particles[h], boxIndex, evaluationTime);
//                        writeElectricFieldToFile(Grid, &particles[0], Fields, step);
//                    }
//                }
//            }
        }
       // *time += dt;
    }
}


void calcFieldsOnGridWithoutNearField(particle *particles, Grid *Grid, Fields *Fields, int numberOfParticles, double t) {
    
//    printf("time at which fields are calculated = %f\n", t);
    
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int p = 0; p < numberOfParticles; p++){
        for (int boxIndex = 0; boxIndex < numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ; boxIndex++){
            if(!boxIsInNearFieldOfParticle(Grid, &particles[p], boxIndex)){
//                addLWFieldsInBox(Grid, &particles[p], boxIndex, t);
                addLWFieldsInBoxforTest(Grid, Fields, &particles[p], boxIndex, t);
            }
        }
    }
}

///@brief In the context of analytial precalculation of fields, a reversed particle pusher is used in order to obtain particle positions x_i, at which particles have to start, in order to be at desired initial conditions after 'precalculationTime'. Given the fact, that particles start moving in 'normal' direction at x_i, variables like iteration count have to be set to 0 again and velocities need to be reversed again.
void resetInitialConditions(particle *particles, int numberOfParticles, double *Bext) {
    scaleVectorByConstantFactor(Bext, -1.0);
    
    for(int p = 0; p < numberOfParticles; p++) {
        for(int i = 0; i < 3; i++) {
            particles[p].uRel[i + 1] *= -1;
            //particles[p].uRelHistory[0][i+1] *= -1;
        }
        particles[p].historyLength = 0;
        particles[p].iterationCount = 0;
        particles[p].xRel[0] = 0;
    }
    //*time = 0;
}

void writeHistoryOfParticlesToFile(particle *particles, char *filename, int actualTimeIndex, int numberOfParticles) {
    for(int p = 0; p < numberOfParticles; p++) {
        FILE *fid = NULL;
        sprintf(filename, "Particle%d_timeStep%d", p, actualTimeIndex);
        strcat(filename, ".txt");
        fid = fopen(filename, "w");
        
        if(fid == NULL) {
            printf("ERROR: Not able to open file\n");
        }
        
        for(int i = 0; i < 4; i++) {
            fprintf(fid, "%f\t", particles[p].xRel[i]);
        }
        fprintf(fid, "\n");
        
        for(int i = 0; i < 4; i++) {
            fprintf(fid, "%f\t", particles[p].uRel[i]);
        }
        fprintf(fid, "\n");
        fclose(fid);
    }
}


///@brief Nyström Algorithm for solving DGL for particle interacting with fields of other particles
///@param p This is the outer for-loop, indicating the actual time iteration
///@param particle Array oof particle structs will be passed in. Within method, actual particle will be accesses by index "actualParticle" within particle-struct
///@param numberOfParticles total number of particles in simulation
void Nystrom(particle *particle, Fields *Fields, Grid *Grid, int p, double dt, int numberOfParticles, int actualParticle, double realdt, int tN, double Bext[3], double Eext[3], bool addHistoryToFile){
    
    double E[3] = {0};
    double B[3] = {0};
    
    // fill E, B with external fields no matter what (even outside of simulation area)
    for(int i = 0; i < 3; i++) {
        E[i] = Eext[i];
        B[i] = Bext[i];
    }
    
    double F[4][4] = {0};
    
    //double sol[n];
    double a[4];
    // double be[4];
    double b[4];
    // double ce[4];
    double c[4];
    // double de[4];
    double d[4];
    double u_tmp[4], x_tmp[4];
   
    // For all Far Fields: Access saved B- and E field at closest grid position (CGP) to particle and pass those Arrays E[3] and B[3] into force tensor together with near field interactions.
    int i = (particle[actualParticle].xRel[1]) / (Grid->dx);
    int j = (particle[actualParticle].xRel[2]) / (Grid->dy);
    int k = (particle[actualParticle].xRel[3]) / (Grid->dz);
    
    int ny = Grid->Ny;
    int nz = Grid->Nz;
    
    // Far fields created by other particles only exist on grid
    if(i >= 0 && j >= 0 && k >= 0) {
        
        E[0] += Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 0];
        E[1] += Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 1];
        E[2] += Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 2];
        
        B[0] += Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 0];
        B[1] += Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 1];
        B[2] += Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 2];
    }
    
    int currentBoxIndexOfOtherParticle = -1;
    for (int d = 0; d < numberOfParticles; d++) {
        if(d != actualParticle) {
            currentBoxIndexOfOtherParticle = calcCurrentBoxIndexOfParticle(&particle[d], Grid);
            if(boxIsInNearFieldOfParticle(Grid, &particle[actualParticle], currentBoxIndexOfOtherParticle)) {
                
                double beta[3] = {0};
                double intersectionPoint[4] = {0};
                double velocityAtIntersectionPoint[4] = {0};
                double gamma_sq;
                double R_sq;
                double R;
                double n[3] = {0};
                double betaDot[3] = {0};
                double dt = 0.5 * Grid->dx;
                double ELW[3] = {0};
                double BLW[3] = {0};

                int currentHistoryLength = particle[d].iterationCount;

                for (int index = 0; index < currentHistoryLength - 1; index ++) {
                    if(isInsideBackwardLightcone(particle[d].xRelHistory[index], particle[actualParticle].xRel) && !isInsideBackwardLightcone(particle[d].xRelHistory[index + 1], particle[actualParticle].xRel)) {
                        calculateIntersectionPointforTest(particle[d].xRelHistory[index], particle[d].xRelHistory[index + 1], particle[d].uRelHistory[index], particle[d].uRelHistory[index + 1], particle[actualParticle].xRel, intersectionPoint, velocityAtIntersectionPoint);
                        calculateBetaforTest(particle[d].xRelHistory[index], particle[d].xRelHistory[index + 1], beta);
                        calculateLienardWiechertParametersforTest(intersectionPoint, particle[actualParticle].xRel, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
                        calculateBetaDotforTest(particle[d].uRelHistory[index], particle[d].uRelHistory[index + 1], dt, betaDot);
                        calcuateLienardWiechertFieldsforTest(gamma_sq, R_sq, R, n, beta, betaDot, 1, ELW, BLW);
                        
                        for (int i = 0; i < 3; i++){
                            E[i] += ELW[i];
                            B[i] += BLW[i];
                        }
                        break;
                    }
                }
            }
        }
    }
    
    F[0][0] = 0;
    F[0][1] = E[0];
    F[0][2] = E[1];
    F[0][3] = E[2];
    F[1][0] = E[0];
    F[1][1] = 0;
    F[1][2] = B[2];
    F[1][3] = -B[1];
    F[2][0] = E[1];
    F[2][1] = -B[2];
    F[2][2] = 0;
    F[2][3] = B[0];
    F[3][0] = E[2];
    F[3][1] = B[1];
    F[3][2] = -B[0];
    F[3][3] = 0;
    
    CalcFForNystrom(a, particle[actualParticle].xRel, particle[actualParticle].uRel, 1, E, B, F);
    
    for(int i=0; i<4; i++){
        x_tmp[i] = particle[actualParticle].xRel[i] + (dt / 2) * particle[actualParticle].uRel[i] + (dt * dt/8) * a[i];
        u_tmp[i] = particle[actualParticle].uRel[i] + dt/2 * a[i];
    }
    
    CalcFForNystrom(b, x_tmp, u_tmp, 1, E, B, F);
    
    for(int i=0; i<4; i++) {
        u_tmp[i] = particle[actualParticle].uRel[i] + (dt/2 * b[i]);
    }
    
    CalcFForNystrom(c, x_tmp, u_tmp, 1, E, B, F);
    
    for (int i=0; i<4; i++){
        x_tmp[i] = particle[actualParticle].xRel[i] + (dt * particle[actualParticle].uRel[i]) + (dt * dt / 2) * c[i];
        u_tmp[i] = particle[actualParticle].uRel[i] + dt * c[i];
    }
    
    CalcFForNystrom(d, x_tmp, u_tmp, 1, E, B, F);
    
    for(int i = 0; i<4; i++) {
        particle[actualParticle].xRel[i] = particle[actualParticle].xRel[i] + dt * particle[actualParticle].uRel[i] + (dt * dt / 6) * (a[i] + b[i] + c[i]);
        particle[actualParticle].uRel[i] = particle[actualParticle].uRel[i] + (dt / 6) * (a[i] + 2 * b[i] + 2 * c[i] + d[i]);
    }
//    particle[actualParticle].xRelHistory[p+1][0] = particle[actualParticle].xRel[0];
//    particle[actualParticle].xRelHistory[p+1][1] = particle[actualParticle].xRel[1];
//    particle[actualParticle].xRelHistory[p+1][2] = particle[actualParticle].xRel[2];
//    particle[actualParticle].xRelHistory[p+1][3] = particle[actualParticle].xRel[3];
//
//    particle[actualParticle].uRelHistory[p+1][0] = particle[actualParticle].uRel[0];
//    particle[actualParticle].uRelHistory[p+1][1] = particle[actualParticle].uRel[1];
//    particle[actualParticle].uRelHistory[p+1][2] = particle[actualParticle].uRel[2];
//    particle[actualParticle].uRelHistory[p+1][3] = particle[actualParticle].uRel[3];
//
//    particle[actualParticle].iterationCount += 1;
}

void addCurrentStateToParticleHistory(particle *particle,  int index) {
    
    for (int i = 0; i < 4; i++){
        particle->xRelHistory[particle->historyLength][i] = particle->xRel[i];
        particle->uRelHistory[particle->historyLength][i] = particle->uRel[i];
    }
    particle->historyLength += 1;
    particle->iterationCount += 1;
}

/// Boris pusher velocity update
void updateVelocityWithBorisPusher(particle *Particles, Grid *Grid, Fields * Fields, int numberOfParticles, int particleIndex, double Eextern[3], double Bextern[3], double dt, int actualTimeStep) {
    
    particle *Particle = &Particles[particleIndex];
    
    int dimension = 3;
    double uPrime[3];
    double t[3];
    double s[3];
    double uMinusCrossT[3];
    double uPrimeCrossS[3];
    double E[3]={0};
    double B[3]={0};
    double absoluteSquareValueOfT;
    double uModifiedForCrossProduct[3];
    double chargeOverMass = Particle->charge / Particle->mass;
 
    // fill E, B with external fields no matter what (even outside of simulation area)
    for(int i = 0; i < 3; i++) {
        E[i] = Eextern[i];
        B[i] = Bextern[i];
    }
    
    // For all Far Fields: Access saved B- and E field at closest grid position (CGP) to particle and pass those Arrays E[3] and B[3] into force tensor together with near field interactions.
    int i = Particle->xRel[1] / (Grid->dx);
    int j = (Particle->xRel[2]) / (Grid->dy);
    int k = (Particle->xRel[3]) / (Grid->dz);
    
    int ny = Grid->Ny;
    int nz = Grid->Nz;
    
    // Far fields created by other particles only exist on grid
    if(i >= 0 && j >= 0 && k >= 0) {
        
        E[0] += Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 0];
        E[1] += Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 1];
        E[2] += Fields->electricField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 2];
        
        B[0] += Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 0];
        B[1] += Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 1];
        B[2] += Fields->magneticField[3 * (nz) * (ny) * (i) + 3 * (nz) * (j) + 3 * (k) + 2];
    }
    
    int currentBoxIndexOfOtherParticle = -1;
    for (int d = 0; d < numberOfParticles; d++) { // ATTENTION: delete this after testing
        if(d != particleIndex) { // && d != 0
            currentBoxIndexOfOtherParticle = calcCurrentBoxIndexOfParticle(&Particles[d], Grid);
            if(boxIsInNearFieldOfParticle(Grid, &Particles[particleIndex], currentBoxIndexOfOtherParticle)) {
                
                double beta[3] = {0};
                double intersectionPoint[4] = {0};
                double velocityAtIntersectionPoint[4] = {0};
                double gamma_sq;
                double R_sq;
                double R;
                double n[3] = {0};
                double betaDot[3] = {0};
                double dt = 0.5 * Grid->dx;
                double ELW[3] = {0};
                double BLW[3] = {0};
                
                int currentHistoryLength = Particles[d].iterationCount;
                
                for (int index = 0; index < currentHistoryLength - 1; index ++) {
                    if(isInsideBackwardLightcone(Particles[d].xRelHistory[index], Particles[particleIndex].xRel) && !isInsideBackwardLightcone(Particles[d].xRelHistory[index + 1], Particles[particleIndex].xRel)) {
                        
                        // get history index, at which fields were emitted (in contrast to actual simulation index)
                        printf("particle history of particle%d\n", d);
                        printf("actual time step = %d\n", actualTimeStep);
                        printf("emitted index = %d\n", index);
                        
                        calculateIntersectionPointforTest(Particles[d].xRelHistory[index], Particles[d].xRelHistory[index + 1], Particles[d].uRelHistory[index], Particles[d].uRelHistory[index + 1], Particles[particleIndex].xRel, intersectionPoint, velocityAtIntersectionPoint);
                        calculateBetaforTest(Particles[d].xRelHistory[index], Particles[d].xRelHistory[index + 1], beta);
                        calculateLienardWiechertParametersforTest(intersectionPoint, Particles[particleIndex].xRel, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
                        calculateBetaDotforTest(Particles[d].uRelHistory[index], Particles[d].uRelHistory[index + 1], dt, betaDot);
                        calcuateLienardWiechertFieldsforTest(gamma_sq, R_sq, R, n, beta, betaDot, 1, ELW, BLW);
                        
                        for (int i = 0; i < 3; i++){
                            E[i] += ELW[i];
                            B[i] += BLW[i];
                        }
                        break;
                    }
                }
            }
        }
    }
    
    // ATTENTION: Delete after testing
    if (particleIndex == 0) {
        // delete all fields (external and interaction fields with other particles), because particle 0 is supposed to be static
        for (int i = 0; i < 3; i++) {
            E[i] = 0;
            B[i] = 0;
        }
    }

    // assisting values
    for (int i = 0; i < dimension; i++){
        t[i] = chargeOverMass * B[i] * dt * 0.5;
    }
    
    absoluteSquareValueOfT = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    
    for(int i = 0; i < dimension; i++){
        s[i] = 2 * t[i] / (1 + absoluteSquareValueOfT);
    }
    
    // implementation of boris method
    // first obtain “uMinus” by adding half acceleration to the initial velocity
    for (int i = 0; i < dimension; i++){
        Particle->uRel[i+1] +=  chargeOverMass * E[i] * dt * 0.5;
    }
    
    // Since u is a four vector rewrite it to a three dimensional vector to make it usabale for crossProduct
    for(int i = 0; i < dimension; i++){
        uModifiedForCrossProduct[i] = Particle->uRel[i+1];
    }
    
    crossProduct(uModifiedForCrossProduct, t, uMinusCrossT);
    
    for (int i = 0; i < dimension; i++){
        uPrime[i] = Particle->uRel[i+1] + uMinusCrossT[i];
    }
    
    crossProduct(uPrime, s, uPrimeCrossS);
    
    // then obtain "uPlus" by performing rotation with "uPrime" and "s"
    for (int i = 0; i < dimension; i++){
        Particle->uRel[i+1] +=  uPrimeCrossS[i];
    }
    
    // finally, add another half acceleration,
    for (int i = 0; i < dimension; i++){
        Particle->uRel[i+1] += chargeOverMass * E[i] * dt * 0.5;
    }
    
    Particle->uRel[0] = getGammaFromVelocityVector(Particle->uRel);
    
    Particles[particleIndex].uRelHistory[actualTimeStep + 1][0] = Particles[particleIndex].uRel[0];
    Particles[particleIndex].uRelHistory[actualTimeStep + 1][1] = Particles[particleIndex].uRel[1];
    Particles[particleIndex].uRelHistory[actualTimeStep + 1][2] = Particles[particleIndex].uRel[2];
    Particles[particleIndex].uRelHistory[actualTimeStep + 1][3] = Particles[particleIndex].uRel[3];
}

void updateLocation(particle *particle, Grid *Grid, double dt, int currentParticle, int currentTimeStep) {
    int dimension = 3;
    
    for(int i = 0; i < dimension; i++) {
        particle[currentParticle].xRel[i+1] += (particle[currentParticle].uRel[i+1] * dt) / particle[currentParticle].uRel[0];
    }
    
    particle[currentParticle].xRel[0] += dt;
    
    particle[currentParticle].xRelHistory[currentTimeStep+1][0] = particle[currentParticle].xRel[0];
    particle[currentParticle].xRelHistory[currentTimeStep+1][1] = particle[currentParticle].xRel[1];
    particle[currentParticle].xRelHistory[currentTimeStep+1][2] = particle[currentParticle].xRel[2];
    particle[currentParticle].xRelHistory[currentTimeStep+1][3] = particle[currentParticle].xRel[3];
    
//    particle[currentParticle].iterationCount += 1;
//    particle[currentParticle].historyLength += 1;
}



///@brief All neccessary operations for one hybrid field push
///@param Grid instance of Grid struct
///@param particle struct containing all particles
///@param numberOfParticles number of particles
///@param t current simulation time
///@param dt time increment
void pushEField(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, double t, double dt) {
    pushEFieldInsideBoxes(Grid, Fields, dt);
    setHFieldOnBorders(Grid,Fields);
    adjustHFields(Grid, Fields, particle, numberOfParticles, t);
    pushEFieldAtBorders(Grid, Fields, dt);
    
}
///@brief All neccessary operations for one hybrid field push
///@param Grid instance of Grid struct
///@param particle struct containing all particles
///@param numberOfParticles number of particles
///@param t current simulation time
///@param dt time increment
void pushHField(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, double t, double dt) {
    pushHFieldInsideBoxes(Grid, Fields, dt);
    setEFieldOnBorders(Grid, Fields);
    adjustEFields(Grid, Fields, particle, numberOfParticles, t);
    pushHFieldAtBorders(Grid, Fields, dt);
}



///@brief Uses Maxwell Pusher to push E-fields inside boxes, except for box borders. This method is used for the hybrid field method.
///@remark E field is calculated using negative curl, i.e. value of H on the left side of E is required. Thus start i,j,k with 1. By default, UPML layer is activated.
///@param Grid Grid struct
///@param dt time increment
void pushEFieldInsideBoxes(Grid *Grid, Fields *Fields, double dt) { //added fields as input parameter
    double Hx_ijk;
    double Hy_ijk;
    double Hz_ijk;
    double Hz_ijm1k;
    double Hy_ijkm1;
    double Hx_ijkm1;
    double Hz_im1jk;
    double Hy_im1jk;
    double Hx_ijm1k;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    double dOld;
    bool useUPML = true;
    
    int numberOfGridPointsInX = Grid->Nx;
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    int i, j, k;
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                
                if (i % Grid->numberOfGridPointsPerBoxInX == 0 || j % Grid->numberOfGridPointsPerBoxInY == 0 || k % Grid->numberOfGridPointsPerBoxInZ == 0){
                    continue;
                } //This is the same code as in pushEField function:
                Hx_ijk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                Hy_ijk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Hz_ijk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                Hz_ijm1k = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j - 1) + 3 * k + 2];
                Hy_ijkm1 = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 1];
                Hx_ijkm1 = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 0];
                Hz_im1jk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                Hy_im1jk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1];
                Hx_ijm1k = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j - 1) + 3 * (k) + 0];
                
                if (useUPML){
                    dOld = Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c1E[j] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c2E[j] * ((Hz_ijk - Hz_ijm1k) / Grid->dy - (Hy_ijk - Hy_ijkm1) / Grid->dz);
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c3E[k] * Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c4E[k] * (Fields->c5E[i] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Fields->c6E[i] * dOld);
                    
                    dOld = Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c1E[k] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c2E[k] * ((Hx_ijk - Hx_ijkm1) / Grid->dz - (Hz_ijk - Hz_im1jk) / Grid->dx); //Factor 0.5
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c3E[i] * Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c4E[i] * (Fields->c5E[j] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Fields->c6E[j] * dOld);
                    
                    
                    dOld = Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c1E[i] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c2E[i] * ((Hy_ijk - Hy_im1jk) / Grid->dx - (Hx_ijk - Hx_ijm1k) / Grid->dy); //Factor 0.5
                    
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c3E[j] * Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c4E[j] * (Fields->c5E[k] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Fields->c6E[k] * dOld);
                    
                } else {
                    double val_x = cny * (Hz_ijk - Hz_ijm1k) - cnz * (Hy_ijk - Hy_ijkm1);
                    double val_y = cnz * (Hx_ijk - Hx_ijkm1) - cnx * (Hz_ijk - Hz_im1jk);
                    double val_z = cnx * (Hy_ijk - Hy_im1jk) - cny * (Hx_ijk - Hx_ijm1k);
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += val_x;
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += val_y;
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += val_z;
                }
            }
        }
    }
}

///@brief for a given "precalculationTime", the nystrom pusher is used to calculate particle history backwards, in order to obtain the particle positions, at which particle has to start, in order to be at defined initial conditions after "precalculationTime". This method is used, if one needs analytic precalculation of fields.
///@param precalculationTime given time span, for which reversed nystrom algorithm is used
///@param externBField vector for external B fields, if given
///@param externEField vector for external E fields, if given
void nystromBackwards(particle *particles, Grid *Grid, Fields *Fields, int numberOfParticles, double dt, double precalculationTime, double *externBField, double *externEField) {
    scaleVectorByConstantFactor(externBField, -1.0);
    
    for(int p = 0; p < numberOfParticles; p++) {
        for(int i = 0; i < 3; i++) {
            particles[p].uRel[i+1] *= -1;
            //particles[p].uRelHistory[0][i+1] *= -1;
        }
    }
    
    for(int preStep = 0; preStep < (precalculationTime / dt); preStep++) {
        for(int p = 0; p < numberOfParticles; p++) {
            addCurrentStateToParticleHistory(&particles[p], preStep);
            
            updateVelocityWithBorisPusher(particles, Grid, Fields, numberOfParticles, p, externEField, externBField, dt, preStep);
            updateLocation(particles, Grid, dt, p, preStep);
            
//            double tau = dt / particles[p].uRel[0];
//            Nystrom(particles, Fields, Grid, preStep, tau, numberOfParticles, p, dt, 20, externBField, externEField, false);
        }
    }
}

///@brief scale (3D) vector by constant factor
void scaleVectorByConstantFactor(double vector[3], double constFactor) {
    for(int i = 0; i < 3; i++) {
        vector[i] *=constFactor;
    }
}

///d@brief in case one wants to precalculate fields analatytically before using the actual hybrid field algorithm, one has to be sure to allocate the additional memory needed for this precalculation
///@param particles test
void reallocateMemoryForParticleHistories(particle *particles, int numberOfParticles, double dt, double simulationTime, double precalculationTime) {
    for(int i = 0; i < numberOfParticles; i++) {
        particles[i].xRelHistory = (double **) realloc(particles[i].xRelHistory, ((simulationTime / dt + precalculationTime / dt) + 1) * sizeof(double*));
        if(particles[i].xRelHistory == NULL) {
            printf("ERROR while trying to reallocate memory for 'xRelHistory' array");
            return;
        }
        for(int j = simulationTime / dt; j <= (simulationTime / dt) + (precalculationTime / dt); j++) {
            particles[i].xRelHistory[j] = (double*) malloc(4 * sizeof(double));
            if(particles[i].xRelHistory[j] == NULL) {
                printf("ERROR while trying to reallocate memory within 'xRelHistory' array in entry %i", j);
                return;
            }
        }
    }
    
    
    for(int i = 0; i < numberOfParticles; i++) {
        particles[i].uRelHistory = (double **) realloc(particles[i].uRelHistory, ((simulationTime / dt + precalculationTime / dt) + 1) * sizeof(double*));
        if(particles[i].uRelHistory == NULL) {
            printf("ERROR while trying to reallocate memory for 'xRelHistory' array");
            return;
        }
        for(int j = simulationTime / dt; j <= (simulationTime / dt) + (precalculationTime / dt); j++) {
            particles[i].uRelHistory[j] = (double*) malloc(4 * sizeof(double));
            if(particles[i].uRelHistory[j] == NULL) {
                printf("ERROR while trying to reallocate memory within 'xRelHistory' array in entry %i", j);
                return;
            }
        }
    }
}



///@brief Uses Maxwell Pusher to push E-fields inside boxes, except for box borders. This method is used for the hybrid field method.
///@remark H field is calculated via positive curl. Therefore value of E on the right side of H on the grid is required. Thus stop i,j,k at n - 1 where n denotes the numberOfGridPoints. UPML is activated by default
///@param Grid instance of Grid struct
///@param dt time increment
void pushHFieldInsideBoxes(Grid *Grid, Fields *Fields, double dt) {
    double Ex_ijk;
    double Ey_ijk;
    double Ez_ijk;
    double Ey_ijkp1;
    double Ez_ijp1k;
    double Ez_ip1jk;
    double Ex_ijkp1;
    double Ex_ijp1k;
    double Ey_ip1jk;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    
    int numberOfGridPointsInX = Grid->Nx;
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    double bOld;
    bool useUPML = true;
    
    int i, j, k;
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                if ((i + 1) % Grid->numberOfGridPointsPerBoxInX == 0 || (j + 1) % Grid->numberOfGridPointsPerBoxInY == 0 || (k + 1) % Grid->numberOfGridPointsPerBoxInZ == 0){
                    continue;
                }
                Ey_ijkp1 = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k + 1) + 1];
                Ey_ijk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1];
                Ez_ijp1k = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j + 1) + 3 * (k) + 2];
                Ez_ijk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2];
                Ez_ip1jk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2];
                Ex_ijkp1 = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k + 1) + 0];
                Ex_ijk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0];
                Ex_ijp1k = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j + 1) + 3 * (k) + 0];
                Ey_ip1jk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1];
                
                if(useUPML){
                    
                    bOld = Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c1H[j] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c2H[j] * ((Ey_ijkp1 - Ey_ijk) / Grid->dz - (Ez_ijp1k - Ez_ijk) / Grid->dy);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c3H[k] * Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c4H[k] * (Fields->c5H[i] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Fields->c6H[i] * bOld);
                    
                    bOld = Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c1H[k] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c2H[k] * ((Ez_ip1jk - Ez_ijk) / Grid->dx - (Ex_ijkp1 - Ex_ijk) / Grid->dz);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c3H[i] * Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c4H[i] * (Fields->c5H[j] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Fields->c6H[j] * bOld);
                    
                    bOld = Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c1H[i] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c2H[i] * ((Ex_ijp1k - Ex_ijk) / Grid->dy - (Ey_ip1jk - Ey_ijk) / Grid->dx);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c3H[j] * Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c4H[j] * (Fields->c5H[k] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Fields->c6H[k] * bOld);
                }
                else{
                    
                    double add_x = cnz * (Ey_ijkp1 - Ey_ijk) - cny * (Ez_ijp1k - Ez_ijk);
                    double add_y = cnx * (Ez_ip1jk - Ez_ijk) - cnz * (Ex_ijkp1 - Ex_ijk);
                    double add_z = cny * (Ex_ijp1k - Ex_ijk) - cnx * (Ey_ip1jk - Ey_ijk);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += add_x;
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += add_y;
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += add_z;
                }
            }
        }
    }
}



///@brief For all boxes, all values of H in the plane left, infront and below of the current box are being stored. Following the Yee-Scheme, only Hy and Hz are needed from the left plane, Hz and Hx are needed from the plane infront and Hy and Hx are needed from the plane below. If the current box is a box where to the left does not exist (ib = 0) then the respective values for H are set to 0.
///@remark Hz_im1 and all others are matrices. The first index denotes the boxIndex and the second Index the gridPoints in the plane.
///@param Grid pointer to Grid struct
void setHFieldOnBorders(Grid *Grid, Fields *Fields) {
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    
    //now use box indicices for looping over all boxes: ib, jb, kb
    for (int ib = 0; ib < numberOfBoxesInX; ib ++){
        for (int jb = 0; jb < numberOfBoxesInY; jb ++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb ++){
                
                if (ib != 0){
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++) {
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Fields->Hz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX - 1) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                            Fields->Hy_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX - 1) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 1];
                        }
                    }
                }
                else{
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Fields->Hz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                            Fields->Hy_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                        }
                    }
                }
                
                if (jb != 0){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Fields->Hx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY - 1) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 0];
                            Fields->Hz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY - 1) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            Fields->Hx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                            Fields->Hz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                        }
                    }
                }
                
                if (kb != 0){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            
                            Fields->Hx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ - 1) + 0];
                            Fields->Hy_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ - 1) + 1];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            Fields->Hx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                            Fields->Hy_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                        }
                    }
                }
            }
            
            
        }
    }
}


///@brief For all boyes all values of E in the plane right, behind and above of the current box are stored. Following the Yee-Algorithm, only Ey and Ez are needed from the right plane. Ez and Ex are needed from the plane behind and Ey and Ex are needed from the plane above. If the current box is a box where the plane to the right does not exist (ib = numberOfBoxesInX) then the respective values for E are set to 0. Used Convention: e.g. Ey_ip1 = y-component of electric field, indicating the grid points in the plane right next to the box in x-direction (i),...etc.
///@remark Ey_ip1 and all others are matrices. The first index denotes the boxIndex (where boxes are counted analogue to grid points) and the second Index the gridPoints in the plane, being indexed by numbers smaller or equal to the number of grid points in the corresponding direction.
///@param Grid pointer to Grid struct
void setEFieldOnBorders(Grid *Grid, Fields *Fields) {
    
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int ib, jb, kb;
    
    for (ib = 0; ib < numberOfBoxesInX; ib++){
        for (jb = 0; jb < numberOfBoxesInY; jb++){
            for (kb = 0; kb < numberOfBoxesInZ; kb++){
                
                if (ib != numberOfBoxesInX - 1){
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Fields->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * ((ib + 1) * numberOfGridPointsForBoxInX) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 1];
                            Fields->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * ((ib + 1) * numberOfGridPointsForBoxInX) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                        }
                    }
                }
                else{
                    for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Fields->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                            Fields->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * jd + kd] = 0;
                        }
                    }
                }
                
                if (jb != numberOfBoxesInY - 1){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Fields->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * ((jb + 1) * numberOfGridPointsForBoxInY) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 2];
                            Fields->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * ((jb + 1) * numberOfGridPointsForBoxInY) + 3 * (kb * numberOfGridPointsForBoxInZ + kd) + 0];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                            
                            Fields->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                            Fields->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * id + kd] = 0;
                        }
                    }
                }
                
                if (kb != numberOfBoxesInZ - 1){
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            
                            Fields->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * ((kb + 1) * numberOfGridPointsForBoxInZ) + 0];
                            Fields->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (ib * numberOfGridPointsForBoxInX + id) + 3 * numberOfGridPointsInZ * (jb * numberOfGridPointsForBoxInY + jd) + 3 * ((kb + 1) * numberOfGridPointsForBoxInZ) + 1];
                        }
                    }
                }
                else{
                    for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                        for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                            
                            Fields->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                            Fields->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * id + jd] = 0;
                        }
                    }
                }
            }
        }
    }
}


///@brief For all boxes, the H fields in the plane to the left, infront and below are being properly adjusted.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param t current simulation time
void adjustHFields(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const double t) {
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int ib = 0; ib < numberOfBoxesInX; ib++){
        for (int jb = 0; jb < numberOfBoxesInY; jb++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb++){
                
                int boxIndex = ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
                
                if (ib != 0){
                    adjustHyz_im1(Grid, Fields, particle, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (jb != 0){
                    adjustHxz_jm1(Grid, Fields, particle, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (kb != 0){
                    adjustHxy_km1(Grid, Fields, particle, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
            }
        }
    }
}


///@brief For all boxes the E fields in the plane to the left, infront and below are being adjusted.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param t current simulation time
void adjustEFields(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const double t) {
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    for (int ib = 0; ib < numberOfBoxesInX; ib++){
        for (int jb = 0; jb < numberOfBoxesInY; jb++){
            for (int kb = 0; kb < numberOfBoxesInZ; kb++){
                
                int boxIndex = ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
                
                if (ib != numberOfBoxesInX - 1){
                    adjustEyz_ip1(Grid, Fields, particle, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (jb != numberOfBoxesInY - 1){
                    adjustExz_jp1(Grid, Fields, particle, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
                
                if (kb != numberOfBoxesInZ - 1){
                    adjustExy_kp1(Grid, Fields, particle, numberOfParticles, boxIndex, ib, jb, kb, t);
                }
            }
        }
    }
}



///@brief Method adjusts the E values on the left and right side border of the near field.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexIp1 is not. Then we are at the right border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value to the right. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the left border of the near field, i.e. current box is the not in near field but boxIndexIp1 is, then we push HField values from the far field. For this we need values to the rigth, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustEyz_ip1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t) {
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int boxIndexIp1 = calcBoxIndexIp1(Grid, boxIndex);
    double xObserver[4];
    
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexIp1)){
    
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * ((ib + 1) * numberOfGridPointsForBoxInX);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    subLWField(Grid, &particle[p], &Fields->Ey_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 1);
                    subLWField(Grid, &particle[p], &Fields->Ez_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 2);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexIp1)){
            
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * ((ib + 1) * numberOfGridPointsForBoxInX);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &particle[p], &(Fields->Ey_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 1);
                    addLWField(Grid, &particle[p], &(Fields->Ez_ip1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 2);
                }
            }
            
        }
    }//end of for Loop for multiple partciles
    
}



///@brief Method adjusts the E values infront and on the back side border of the near field.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexJp1 is not. Then we are at the back side border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value behind. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the front side border of the near field, i.e. current box is the not in near field but boxIndexJp1 is, then we push HField values from the far field. For this we need values infront, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustExz_jp1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t) {
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int boxIndexJp1 = calcBoxIndexJp1(Grid, boxIndex);
    double xObserver[4];
    
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexJp1)){
    
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * ((jb + 1) * numberOfGridPointsForBoxInY);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    subLWField(Grid, &particle[p], &Fields->Ex_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 0);
                    subLWField(Grid, &particle[p], &Fields->Ez_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 2);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexJp1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * ((jb + 1) * numberOfGridPointsForBoxInY);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &particle[p], &(Fields->Ex_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 0);
                    addLWField(Grid, &particle[p], &(Fields->Ez_jp1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 2);
                }
            }
            
        }
    }
}


///@brief Method adjusts the E values at top and bottom side border of the near field.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustEFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexKp1 is not. Then we are at the top side border of the near field. Since we want to push a value for the Hfield inside the near field we need (among others) the E field value above. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the bottom side border of the near field, i.e. current box is the not in near field but boxIndexKp1 is, then we push HField values from the far field. For this we need values above, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the H field value in the far field correctly.
void adjustExy_kp1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t) {
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int boxIndexKp1 = calcBoxIndexKp1(Grid, boxIndex);
    double xObserver[4];
    
    for(int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexKp1)){
    
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * ((kb + 1) * numberOfGridPointsForBoxInZ);
                    
                    subLWField(Grid, &particle[p], &Fields->Ex_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 0);
                    subLWField(Grid, &particle[p], &Fields->Ey_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 1);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexKp1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * ((kb + 1) * numberOfGridPointsForBoxInZ);
                    
                    addLWField(Grid, &particle[p], &(Fields->Ex_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 0);
                    addLWField(Grid, &particle[p], &(Fields->Ey_kp1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 1);
                }
            }
            
        }
    }
}


///@brief Method adjusts the H values on the left and right side border of the near field.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexIm1 is not. Then we are at the left border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value to the left. Since this point is in the far field, the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the right border of the near field, i.e. current box is the not in near field but boxIndexIm1 is, then we push EField values from the far field. For this we need values to the left, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHyz_im1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t) {
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    //Calculate the box index to the left in x direction to the current box index
    int boxIndexIm1 = calcBoxIndexIm1(Grid, boxIndex);
    //xobserver, which is looping through the whole surface on which fields need potentially to be adjusted
    double xObserver[4];
    
    // This is for multiple particles
    for (int p = 0; p < numberOfParticles; p++) {
        if (boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexIm1)){
    
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t; //absolute time, not index --> was deleted, cause not needed
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX - 1);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    subLWField(Grid, &particle[p], &Fields->Hy_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 4);
                    subLWField(Grid, &particle[p], &Fields->Hz_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd], xObserver, 5);
                }
            }
        }
        else if (!boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexIm1)){
            
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX - 1);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &particle[p], &(Fields->Hy_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 4);
                    addLWField(Grid, &particle[p], &(Fields->Hz_im1[boxIndex][numberOfGridPointsForBoxInZ * jd + kd]), xObserver, 5);
                }
            }
            
        }
    }
}


///@brief Method adjusts the H values infront and on the back side border of the near field.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexJm1 is not. Then we are at the facing border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value infront. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the back border of the near field, i.e. current box is the not in near field but boxIndexJm1 is, then we push EField values from the far field. For this we need values infront, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHxz_jm1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t) {
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int boxIndexJm1 = calcBoxIndexJm1(Grid, boxIndex);
    double xObserver[4];
    
    for (int p = 0; p < numberOfParticles; p++){
            if (boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
                && !boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexJm1)){
    
        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY - 1);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                
                subLWField(Grid, &particle[p], &Fields->Hx_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 3);
                subLWField(Grid, &particle[p], &Fields->Hz_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd], xObserver, 5);
            }
        }
    }
        else if (!boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexJm1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int kd = 0; kd < numberOfGridPointsForBoxInZ; kd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY - 1);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ + kd);
                    
                    addLWField(Grid, &particle[p], &(Fields->Hx_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 3);
                    addLWField(Grid, &particle[p], &(Fields->Hz_jm1[boxIndex][numberOfGridPointsForBoxInZ * id + kd]), xObserver, 5);
                }
            }
            
        }
    }
}



///@brief Method adjusts the H values on the top and bottom side border of the near field.
///@param Grid pointer to Grid struct
///@param particle pointer to Particles struct array
///@param numberOfParticles number of particles
///@param boxIndex current boxIndex from outer loop
///@param ib boxIndex in x direction
///@param jb boxIndex in y direction
///@param kb boxIndex in z direction
///@param t current simulation time
///@remark this method gets called over and over again from "adjustHFields()" method while looping through all boxes. We  check if the current box with boxIndex is in the near field of a particle and the box with boxIndexKm1 is not. Then we are at bottom border of the near field. Since we want to push a value for the Efield inside the near field we need (among others) the value below. Since this point is in the far field the LW fields are already stored on that grid point. Therefore we need to substract the LW field to push the value inside the near field correctly. On the other hand, if we are on the top border of the near field, i.e. current box is the not in near field but boxIndexKm1 is, then we push EField values from the far field. For this we need values below, i.e. in the near field. Since in the near field area no LW fields are stored, we need to calculate them and add them to push the E field value in the far field correctly.
void adjustHxy_km1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t) {
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int boxIndexKm1 = calcBoxIndexKm1(Grid, boxIndex);
    double xObserver[4];
    
    for (int p = 0; p < numberOfParticles; p++){
        if (boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
            && !boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexKm1)){

        for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
            for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                
                xObserver[0] = t;
                xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ - 1);
                
                subLWField(Grid, &particle[p], &Fields->Hx_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 3);
                subLWField(Grid, &particle[p], &Fields->Hy_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd], xObserver, 4);
            }
        }
    }
        else if (!boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndex)
                 && boxIsInNearFieldOfParticle(Grid, &particle[p], boxIndexKm1)){
            
            for (int id = 0; id < numberOfGridPointsForBoxInX; id++){
                for (int jd = 0; jd < numberOfGridPointsForBoxInY; jd++){
                    
                    xObserver[0] = t;
                    xObserver[1] = (Grid->dx) * (ib * numberOfGridPointsForBoxInX + id);
                    xObserver[2] = (Grid->dy) * (jb * numberOfGridPointsForBoxInY + jd);
                    xObserver[3] = (Grid->dz) * (kb * numberOfGridPointsForBoxInZ - 1);
                    
                    addLWField(Grid, &particle[p], &(Fields->Hx_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 3);
                    addLWField(Grid, &particle[p], &(Fields->Hy_km1[boxIndex][numberOfGridPointsForBoxInY * id + jd]), xObserver, 4);
                }
            }
            
        }
    }
}



///@brief Pushes EField at box borders using the Maxwell-Pusher.
///@remark E field is calculated via negative curl. Therefore value of H on the left side of E on the grid is required. Thus start i,j,k with 1. UPML is activated by default
///@param Grid pointer to Grid struct
///@param dt time increment
void pushEFieldAtBorders(Grid *Grid, Fields *Fields, double dt) {
    
    double Hx_ijk;
    double Hy_ijk;
    double Hz_ijk;
    double Hz_ijm1k;
    double Hy_ijkm1;
    double Hx_ijkm1;
    double Hz_im1jk;
    double Hy_im1jk;
    double Hx_ijm1k;
    
    double dOld;
    bool useUPML = true;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    int numberOfGridPointsInX = Grid->Nx;
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int i, j, k;
    
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                
                if (i % numberOfGridPointsForBoxInX != 0 && j % numberOfGridPointsForBoxInY != 0 && k % numberOfGridPointsForBoxInZ != 0){
                    continue;
                }
                
                int ib = i / numberOfGridPointsForBoxInX;
                int jb = j / numberOfGridPointsForBoxInY;
                int kb = k / numberOfGridPointsForBoxInZ;
                
                Hx_ijk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                Hy_ijk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Hz_ijk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                
                Hx_ijm1k = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j - 1) + 3 * k + 0];
                Hy_im1jk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Hz_im1jk = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i - 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                Hx_ijkm1 = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 0];
                Hy_ijkm1 = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k - 1) + 1];
                Hz_ijm1k = Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j - 1) + 3 * k + 2];
                
                
                /*adjust values*/
                if (i % numberOfGridPointsForBoxInX == 0){
                    Hz_im1jk = Fields->Hz_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                    Hy_im1jk = Fields->Hy_im1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                }
                
                if (j % numberOfGridPointsForBoxInY == 0){
                    Hx_ijm1k = Fields->Hx_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                    Hz_ijm1k = Fields->Hz_jm1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                }
                
                if (k % numberOfGridPointsForBoxInZ == 0){
                    Hx_ijkm1 = Fields->Hx_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                    Hy_ijkm1 = Fields->Hy_km1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                }
                
                if (useUPML){
                    dOld = Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c1E[j] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c2E[j] * ((Hz_ijk - Hz_ijm1k) / Grid->dy - (Hy_ijk - Hy_ijkm1) / Grid->dz);
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c3E[k] * Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c4E[k] * (Fields->c5E[i] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Fields->c6E[i] * dOld);
                    
                    dOld = Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c1E[k] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c2E[k] * ((Hx_ijk - Hx_ijkm1) / Grid->dz - (Hz_ijk - Hz_im1jk) / Grid->dx); //Factor 0.5
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c3E[i] * Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c4E[i] * (Fields->c5E[j] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Fields->c6E[j] * dOld);
                    
                    dOld = Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c1E[i] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c2E[i] * ((Hy_ijk - Hy_im1jk) / Grid->dx - (Hx_ijk - Hx_ijm1k) / Grid->dy); //Factor 0.5
                    
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c3E[j] * Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c4E[j] * (Fields->c5E[k] * Fields->dField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Fields->c6E[k] * dOld);
                    
                }
                else{
                    double val_x = cny * (Hz_ijk - Hz_ijm1k) - cnz * (Hy_ijk - Hy_ijkm1);
                    double val_y = cnz * (Hx_ijk - Hx_ijkm1) - cnx * (Hz_ijk - Hz_im1jk);
                    double val_z = cnx * (Hy_ijk - Hy_im1jk) - cny * (Hx_ijk - Hx_ijm1k);
                    
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += val_x;
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += val_y;
                    Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += val_z;
                }
            }
        }
    }
}


///@brief Pushes HField at box borders using the Maxwell-Pusher.
///@remark H field is calculated via positive curl. Therefore value of E on the right side of H on the grid is required. Thus stop i,j,k at n - 1 where n denotes the numberOfGridPoints. UPML is activated by default
///@param Grid pointer to Grid struct
///@param dt time increment
void pushHFieldAtBorders(Grid *Grid, Fields *Fields, double dt) {
    
    double Ex_ijk;
    double Ey_ijk;
    double Ez_ijk;
    double Ey_ijkp1;
    double Ez_ijp1k;
    double Ez_ip1jk;
    double Ex_ijkp1;
    double Ex_ijp1k;
    double Ey_ip1jk;
    
    //Only used for UMPL
    double bOld;
    bool useUPML = true;
    
    double cnx = 0.5 * dt / Grid->dx;
    double cny = 0.5 * dt / Grid->dy;
    double cnz = 0.5 * dt / Grid->dz;
    
    int numberOfGridPointsInX = Grid->Nx;
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    int i, j, k;
    
    for (i = 1; i < numberOfGridPointsInX - 1; i++){
        for (j = 1; j < numberOfGridPointsInY - 1; j++){
            for (k = 1; k < numberOfGridPointsInZ - 1; k++){
                
                if ((i + 1) % numberOfGridPointsForBoxInX != 0 && (j + 1) % numberOfGridPointsForBoxInY != 0 && (k + 1) % numberOfGridPointsForBoxInZ != 0){
                    continue;
                }
                
                int ib = i / numberOfGridPointsForBoxInX;
                int jb = j / numberOfGridPointsForBoxInY;
                int kb = k / numberOfGridPointsForBoxInZ;
                
                Ex_ijk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                Ey_ijk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Ez_ijk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                
                Ez_ip1jk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                Ey_ip1jk = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i + 1) + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                Ex_ijp1k = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j + 1) + 3 * k + 0];
                Ez_ijp1k = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * (j + 1) + 3 * k + 2];
                Ex_ijkp1 = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k + 1) + 0];
                Ey_ijkp1 = Fields->electricField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * (k + 1) + 1];
                
                /*adjust values*/
                if ((i + 1) % numberOfGridPointsForBoxInX == 0){
                    Ez_ip1jk =
                    Fields->Ez_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                    Ey_ip1jk =
                    Fields->Ey_ip1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (j % numberOfGridPointsForBoxInY) + (k % numberOfGridPointsForBoxInZ)];
                    
                }
                
                if ((j + 1) % numberOfGridPointsForBoxInY == 0){
                    Ex_ijp1k =
                    Fields->Ex_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                    Ez_ijp1k =
                    Fields->Ez_jp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInZ * (i % numberOfGridPointsForBoxInX) + (k % numberOfGridPointsForBoxInZ)];
                }
                
                if ((k + 1) % numberOfGridPointsForBoxInZ == 0){
                    Ex_ijkp1 =
                    Fields->Ex_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                    Ey_ijkp1 =
                    Fields->Ey_kp1[numberOfBoxesInZ * numberOfBoxesInY * ib + numberOfBoxesInZ * jb + kb][numberOfGridPointsForBoxInY * (i % numberOfGridPointsForBoxInX) + (j % numberOfGridPointsForBoxInY)];
                }
                
                if(useUPML){
                    
                    bOld = Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0];
                    
                    Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c1H[j] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c2H[j] * ((Ey_ijkp1 - Ey_ijk) / Grid->dz - (Ez_ijp1k - Ez_ijk) / Grid->dy);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] = Fields->c3H[k] * Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] + Fields->c4H[k] * (Fields->c5H[i] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 0] - Fields->c6H[i] * bOld);
                    
                    bOld = Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1];
                    
                    Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c1H[k] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c2H[k] * ((Ez_ip1jk - Ez_ijk) / Grid->dx - (Ex_ijkp1 - Ex_ijk) / Grid->dz);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] = Fields->c3H[i] * Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] + Fields->c4H[i] * (Fields->c5H[j] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 1] - Fields->c6H[j] * bOld);
                    
                    bOld = Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2];
                    
                    Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c1H[i] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c2H[i] * ((Ex_ijp1k - Ex_ijk) / Grid->dy - (Ey_ip1jk - Ey_ijk) / Grid->dx);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] = Fields->c3H[j] * Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] + Fields->c4H[j] * (Fields->c5H[k] * Fields->bField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * i + 3 * numberOfGridPointsInZ * j + 3 * k + 2] - Fields->c6H[k] * bOld);
                }
                else{
                    
                    double val_x = cnz * (Ey_ijkp1 - Ey_ijk) - cny * (Ez_ijp1k - Ez_ijk);
                    double val_y = cnx * (Ez_ip1jk - Ez_ijk) - cnz * (Ex_ijkp1 - Ex_ijk);
                    double val_z = cny * (Ex_ijp1k - Ex_ijk) - cnx * (Ey_ip1jk - Ey_ijk);
                    
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 0] += val_x;
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 1] += val_y;
                    Fields->magneticField[3 * numberOfGridPointsInZ * numberOfGridPointsInY * (i) + 3 * numberOfGridPointsInZ * (j) + 3 * (k) + 2] += val_z;
                }
            }
        }
    }
}


//===============================================================================================
// Helper Functions for Hybrid Fields
//===============================================================================================

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex to the left of the current box where particle is located.
int calcBoxIndexIm1(Grid *Grid, const int boxIndex) {
    return boxIndex - (Grid->numberOfBoxesInZ) * (Grid->numberOfBoxesInY);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex in front of the current box where particle is located.
int calcBoxIndexJm1(Grid *Grid, const int boxIndex) {
    return boxIndex - (Grid->numberOfBoxesInZ);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex below of the current box where particle is located.
int calcBoxIndexKm1(Grid *Grid, const int boxIndex) {
    return boxIndex - 1;
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex to the right of the current box where particle is located.
int calcBoxIndexIp1(Grid *Grid, const int boxIndex) {
    return boxIndex + (Grid->numberOfBoxesInZ) * (Grid->numberOfBoxesInY);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex behind the current box where particle is located.
int calcBoxIndexJp1(Grid *Grid, const int boxIndex) {
    return boxIndex + (Grid->numberOfBoxesInZ);
}

///@param Grid pointer to Grid struct.
///@param boxIndex current boxIndex of particle
///@returns the boxIndex above the current box where particle is located.
int calcBoxIndexKp1(Grid *Grid, const int boxIndex) {
    return boxIndex + 1;
}


///@brief Bool, if "boxindex" is within Near-Field of particle. Near Fiels is defined as all adjacent boxes to box, in which particle is locared.
///@param Grid pointer to an instance of a Grid struct.
///@param particle pointer to an instance of a Particle struct
///@param boxIndex current box index from outter loop
///@return true if box with Index "boxIndex" is in near field of particle or false if not.
bool boxIsInNearFieldOfParticle(Grid *Grid, particle *particle, int boxIndex){
    int boxIndizesOfNextNeighbourBoxes[27] = {0};
    calcBoxIndizesOfNextNeighbourBoxes(Grid, particle, boxIndizesOfNextNeighbourBoxes);
    
    bool isInNF = false;
    
    for (int i = 0; i < 27; i++){
        if(boxIndizesOfNextNeighbourBoxes[i] == boxIndex){
            isInNF = true;
            break;
        }
    }
    return isInNF;
}



///@brief Calculated the box indices of all boxes adjacent to box, in which particle is located
///@param Grid pointer to an instance of a Grid struct.
///@param particle pointer to an instance of a Particle struct
///@param boxIndizesOfNextNeighbourBoxes array in which the indicies will be stored.
void calcBoxIndizesOfNextNeighbourBoxes(Grid *Grid, particle *particle, int boxIndizesOfNextNeighbourBoxes[27]) {
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    int numberOfBoxesInX = Grid->numberOfBoxesInX;
    
    int currentBoxIndexOfParticle = calcCurrentBoxIndexOfParticle(particle, Grid);
    int indexOfNextNeighbourBox = 0;
    int maxBoxIndex = numberOfBoxesInX * numberOfBoxesInY * numberOfBoxesInZ;
    int index = 0;
    
    for (int ii = -1; ii <= 1; ii++){
        for(int ij = -1; ij <= 1; ij++){
            for(int ik = -1; ik <= 1; ik++){
                indexOfNextNeighbourBox = currentBoxIndexOfParticle + (ik + numberOfBoxesInZ * ij + numberOfBoxesInZ * numberOfBoxesInY * ii);
                if (indexOfNextNeighbourBox < 0 || indexOfNextNeighbourBox > maxBoxIndex){
                    boxIndizesOfNextNeighbourBoxes[index] = -1;
                }
                else{
                    boxIndizesOfNextNeighbourBoxes[index] = indexOfNextNeighbourBox;
                }
                index++;
            }
        }
        
    }
}


///@brief calcualtes the box index in which the particle is currently located.
///@param Grid pointer to an instance of a Grid struct.
///@param particle pointer to an instance of a Particle struct
///@return index of box where particle is currenty located.
///@remark index for box is counted the same way as for E and B field array. First z component, then y and then x.
int calcCurrentBoxIndexOfParticle(particle *particle, Grid *Grid) {
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->dx;
    double dy = Grid->dy;
    double dz = Grid->dz;
    
    int ib = particle->xRel[1] / dx / numberOfGridPointsForBoxInX;
    int jb = particle->xRel[2] / dy / numberOfGridPointsForBoxInY;
    int kb = particle->xRel[3] / dz / numberOfGridPointsForBoxInZ;
    
    return ib * numberOfBoxesInY * numberOfBoxesInZ + jb * numberOfBoxesInZ + kb;
}


///@brief At a certain observation point "xObserver" (will be position of particle), LW fields are being calculated and then added to "destination", i.e. already given field value.
///@param Grid pointer to an instance of a Grid struct.
///@param xObserver (absolute) observation point at which the LW fields shall be calculted, i.e. indices i,j,k need to be multiplicated by respective grid resolution dx, dy, dz.
///@param component component in which the observation point shall be shifted, where 0 = E_x, 1 = E_y, 2 = E_z, 3 = H_x, 4 = H_y, 5 = H_z
///@param particle pointer to an instance of a Particle struct
///@param destination pointer to where LW fields shall be added, i.e. the corresponding field entry (E or H) for a observation point i,j,k. For a given observation point: E[i,j,k] = destination
///@remark The Yee-scheme requires that the E and B fields are shifted against each other, which is being taken care of by creating an "xObserverCopy" Array within this methoid. Attention: Lee-ALgorithm only holds for field calculatiions.
void addLWField(Grid *Grid, particle *particle, double *destination, double xObserver[4], int component) {
    double xObserverShifted[4];
    memcpy(xObserverShifted, xObserver, 4 * sizeof(double));
    switch (component)
    {
        case 0:
            xObserverShifted[1] += 0.5 * (Grid->dx);
            break;
        case 1:
            xObserverShifted[2] += 0.5 * (Grid->dy);
            break;
        case 2:
            xObserverShifted[3] += 0.5 * (Grid->dz);
            break;
        case 3:
            xObserverShifted[2] += 0.5 * (Grid->dy);
            xObserverShifted[3] += 0.5 * (Grid->dz);
            break;
        case 4:
            xObserverShifted[3] += 0.5 * (Grid->dz);
            xObserverShifted[1] += 0.5 * (Grid->dx);
            break;
        case 5:
            xObserverShifted[1] += 0.5 * (Grid->dx);
            xObserverShifted[2] += 0.5 * (Grid->dy);
            break;
    }
    
    double beta[3] = {0};
    double intersectionPoint[4] = {0};
    double velocityAtIntersectionPoint[4] = {0};
    double gamma_sq;
    double R_sq;
    double R;
    double n[3] = {0};
    double betaDot[3] = {0};
    double dt = 0.5 * Grid->dx;
    double E[3] = {0};
    double B[3] = {0};
    
    int currentHistoryLength = particle->iterationCount;
    
    for (int index = 0; index < currentHistoryLength - 1; index ++){
        if(isInsideBackwardLightcone(particle->xRelHistory[index], xObserverShifted) && !isInsideBackwardLightcone(particle->xRelHistory[index + 1], xObserverShifted)) {
            calculateIntersectionPointforTest(particle->xRelHistory[index], particle->xRelHistory[index + 1], particle->uRelHistory[index], particle->uRelHistory[index + 1], xObserverShifted, intersectionPoint, velocityAtIntersectionPoint);
            calculateBetaforTest(particle->xRelHistory[index], particle->xRelHistory[index + 1], beta);
            calculateLienardWiechertParametersforTest(intersectionPoint, xObserverShifted, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
            calculateBetaDotforTest(particle->uRelHistory[index], particle->uRelHistory[index + 1], dt, betaDot);
            calcuateLienardWiechertFieldsforTest(gamma_sq, R_sq, R, n, beta, betaDot, 1, E, B);
            
//            printf("E[0] = %f\n", E[0]);
//            printf("E[1] = %f\n", E[1]);
//            printf("E[2] = %f\n", E[2]);
//            printf("Grid point x = %f\n", xObserver[1]);
//            printf("Grid point y = %f\n", xObserver[2]);
//            printf("Grid point z = %f\n", xObserver[2]);
            
            break;
        }
    }
    
    switch (component)
    {
        case 0:
            *destination += E[0];
            break;
        case 1:
            *destination += E[1];
            break;
        case 2:
            *destination += E[2];
            break;
        case 3:
            *destination += B[0];
            break;
        case 4:
            *destination += B[1];
            break;
        case 5:
            *destination += B[2];
            break;
    }
}


///@brief At a certain observation point "xObserver" (will be position of particle), LW fields are being calculated and then added to "destination", i.e. already given field value.
///@remark This is exactly the same function as "addLWFields" apart from the last part, where field get substracted instead of added
///@param Grid pointer to an instance of a Grid struct.
///@param xObserver (absolute) observation point at which the LW fields shall be calculted, i.e. indices i,j,k need to be multiplicated by respective grid resolution dx, dy, dz.
///@param component component in which the observation point shall be shifted, where 0 = E_x, 1 = E_y, 2 = E_z, 3 = H_x, 4 = H_y, 5 = H_z
///@param particle pointer to an instance of a Particle struct
///@param destination pointer to where LW fields shall ne added.
///@remark The Yee-scheme requires that the E and B fields are shifted against each other, which is being taken care of by creating an "xObserverCopy" Array within this methoid. Attention: Lee-ALgorithm only holds for field calculatiions.
void subLWField(Grid *Grid, particle *particle, double *destination, double xObserver[4], int component) {
    double xObserverShifted[4];
    memcpy(xObserverShifted, xObserver, 4 * sizeof(double));
    switch (component)
    {
        case 0:
            xObserverShifted[1] += 0.5 * (Grid->dx);
            break;
        case 1:
            xObserverShifted[2] += 0.5 * (Grid->dy);
            break;
        case 2:
            xObserverShifted[3] += 0.5 * (Grid->dz);
            break;
        case 3:
            xObserverShifted[2] += 0.5 * (Grid->dy);
            xObserverShifted[3] += 0.5 * (Grid->dz);
            break;
        case 4:
            xObserverShifted[3] += 0.5 * (Grid->dz);
            xObserverShifted[1] += 0.5 * (Grid->dx);
            break;
        case 5:
            xObserverShifted[1] += 0.5 * (Grid->dx);
            xObserverShifted[2] += 0.5 * (Grid->dy);
            break;
    }
    
    double beta[3] = {0};
    double intersectionPoint[4] = {0};
    double velocityAtIntersectionPoint[4] = {0};
    double gamma_sq;
    double R_sq;
    double R;
    double n[3] = {0};
    double betaDot[3] = {0};
    double dt = 0.5 * Grid->dx;
    double E[3] = {0};
    double B[3] = {0};
    
    int currentHistoryLength = particle->iterationCount;
    
    for (int index = 0; index < currentHistoryLength - 1; index ++){
        if(isInsideBackwardLightcone(particle->xRelHistory[index], xObserverShifted) && !isInsideBackwardLightcone(particle->xRelHistory[index + 1], xObserverShifted)) {
            calculateIntersectionPointforTest(particle->xRelHistory[index], particle->xRelHistory[index + 1], particle->uRelHistory[index], particle->uRelHistory[index + 1], xObserverShifted, intersectionPoint, velocityAtIntersectionPoint);
            calculateBetaforTest(particle->xRelHistory[index], particle->xRelHistory[index + 1], beta);
            calculateLienardWiechertParametersforTest(intersectionPoint, xObserverShifted, velocityAtIntersectionPoint, &gamma_sq, &R_sq, &R, n, beta);
            calculateBetaDotforTest(particle->uRelHistory[index], particle->uRelHistory[index + 1], dt, betaDot);
            calcuateLienardWiechertFieldsforTest(gamma_sq, R_sq, R, n, beta, betaDot, 1, E, B);
            break;
        }
    }
    
    switch (component)
    {
        case 0:
            *destination -= E[0];
            break;
        case 1:
            *destination -= E[1];
            break;
        case 2:
            *destination -= E[2];
            break;
        case 3:
            *destination -= B[0];
            break;
        case 4:
            *destination -= B[1];
            break;
        case 5:
            *destination -= B[2];
            break;
    }
}

///@brief Cofficients for UPML Layer are calcuated with a fixed layer width (10 works pretty well).
void precalculateCoefficientsForUPML(Grid *Grid, Fields *Fields) {
    //c1E,H and c2E,H (having y components of sigma and k) as well as c1H and c2H
    for(int j=0;j<Grid->Ny;j++) {
        // How many grid points does the PMP layer contain?
        double layerWidth = 10;
        double dt = 0.25 * Grid->dx;
        double m = 3.5;
        double kmax = 1;
        double sigmaMax = 18; //optimum: \sigma_{max} = 0,8*(m+1)/(\Delta), with \Delta = x,y,z spacing
        double sigmayE;
        double kyE;
        double sigmayH;
        double kyH;
        //left layer
        if(j < layerWidth) {
            //c1,2E: no shift
            sigmayE = pow(((layerWidth-j) / layerWidth),m) * sigmaMax;
            //printf("sigmayE[%d]=%f\n",j,sigmayE);
            kyE = 1+(kmax-1)*pow(((layerWidth-j)/layerWidth),m);
            //c1,2H: y shift of 0,5
            sigmayH = pow(((layerWidth - j - 0.5) / layerWidth),m) * sigmaMax;
            kyH = 1+(kmax-1)*pow(((layerWidth-j-0.5)/layerWidth),m);
            //c1,2E
            Fields->c1E[j] = (2*kyE-sigmayE*dt)/(2*kyE+sigmayE*dt);
            Fields->c2E[j] = 2*dt/(2*kyE+sigmayE*dt);
            //c1,2H
            Fields->c1H[j] = (2*kyH-sigmayH*dt)/(2*kyH+sigmayH*dt);
            Fields->c2H[j] = 2*dt/(2*kyH+sigmayH*dt);
        }
        // simulation area
        if(j >= layerWidth && j < Grid->Ny-layerWidth) {
            sigmayE = 0;
            kyE = 1;
            sigmayH = 0;
            kyH = 1;
            Fields->c1E[j] = (2*kyE-sigmayE*dt)/(2*kyE+sigmayE*dt);
            Fields->c2E[j] = 2*dt/(2*kyE+sigmayE*dt);
            Fields->c1H[j] = (2*kyH-sigmayH*dt)/(2*kyH+sigmayH*dt);
            Fields->c2H[j] = 2*dt/(2*kyH+sigmayH*dt);
        }
        // right layer
        if((Grid->Ny-layerWidth-1) <=j && j < Grid->Ny) {
            sigmayE = pow(((j-Grid->Ny+layerWidth+1)/layerWidth),m)*sigmaMax;
            kyE = 1+(kmax-1)*pow(((j-(Grid->Ny-layerWidth)+1)/layerWidth),m);
            sigmayH = pow(((j+0.5-(Grid->Ny-layerWidth)+1)/layerWidth),m)*sigmaMax;
            kyH = 1+(kmax-1)*pow(((j+0.5-(Grid->Ny-layerWidth)+1)/layerWidth),m);
            Fields->c1E[j] = (2*kyE-sigmayE*dt)/(2*kyE+sigmayE*dt);
            Fields->c2E[j] = 2*dt/(2*kyE+sigmayE*dt);
            Fields->c1H[j] = (2*kyH-sigmayH*dt)/(2*kyH+sigmayH*dt);
            Fields->c2H[j] = 2*dt/(2*kyH+sigmayH*dt);
        }
    }
    
    //c3E,H and c4E,H
    for(int k=0;k<Grid->Nz;k++) {
        double layerWidth = 10;
        double m = 3.5;
        double dt = 0.25*Grid->dx;
        double kmax = 1;
        double sigmaMax = 18;
        double sigmazE;
        double kzE;
        double sigmazH;
        double kzH;
        double epsilon = 1; //Does that need to be changed?
        //left layer
        if(k<layerWidth) {
            //parameters for c3,4E: no z shift
            sigmazE = pow(((layerWidth-k)/layerWidth),m)*sigmaMax;
            //printf("sigmazE[%d]=%f\n",k,sigmazE);
            kzE = 1+(kmax-1)*pow(((layerWidth-k)/layerWidth),m);
            //parameters for c3,4H: z shift of 0,5
            sigmazH = pow(((layerWidth-k-0.5)/layerWidth),m)*sigmaMax;
            kzH = 1+(kmax-1)*pow(((layerWidth-k-0.5)/layerWidth),m);
            //c3,4E
            Fields->c3E[k] = (2*kzE-sigmazE*dt)/(2*kzE+sigmazE*dt);
            Fields->c4E[k] = 1/((2*kzE+sigmazE*dt)*epsilon);
            //c3,4H
            Fields->c3H[k] = (2*kzH-sigmazH*dt)/(2*kzH+sigmazH*dt);
            Fields->c4H[k] = 1/((2*kzH+sigmazH*dt)*epsilon);
        }
        //simulation area
        if(layerWidth<=k && k<(Grid->Nz)-layerWidth) {
            sigmazE = 0;
            kzE = 1;
            sigmazH = 0;
            kzH = 1;
            Fields->c3E[k] = (2*kzE-sigmazE*dt)/(2*kzE+sigmazE*dt);
            Fields->c4E[k] = 1/((2*kzE+sigmazE*dt)*epsilon);
            Fields->c3H[k] = (2*kzH-sigmazH*dt)/(2*kzH+sigmazH*dt);
            Fields->c4H[k] = 1/((2*kzH+sigmazH*dt)*epsilon);
        }
        //right layer
        if((Grid->Nz-layerWidth-1)<=k && k<Grid->Nz) {
            sigmazE = pow(((k-(Grid->Ny-layerWidth)+1)/layerWidth),m)*sigmaMax;
            kzE = 1+(kmax-1)*pow(((k-(Grid->Ny-layerWidth)+1)/layerWidth),m);
            sigmazH = pow(((k+0.5-(Grid->Ny-layerWidth)+1)/layerWidth),m)*sigmaMax;
            kzH = 1+(kmax-1)*pow(((k+0.5-(Grid->Ny-layerWidth)+1)/layerWidth),m);
            Fields->c3E[k] = (2*kzE-sigmazE*dt)/(2*kzE+sigmazE*dt);
            Fields->c4E[k] = 1/((2*kzE+sigmazE*dt)*epsilon);
            Fields->c3H[k] = (2*kzH-sigmazH*dt)/(2*kzH+sigmazH*dt);
            Fields->c4H[k] = 1/((2*kzH+sigmazH*dt)*epsilon);
        }
    }
    
    //c5E,H and c6E,H
    for(int i=0;i<Grid->Nx;i++) {
        double layerWidth = 10;
        double m = 3.5;
        double dt = 0.25*Grid->dx;
        double kmax = 1;
        double sigmaMax = 18;
        double sigmaxE;
        double kxE;
        double sigmaxH;
        double kxH;
        //left layer
        if(i<layerWidth) {
            //parameters for c5,6E: x shift of 0.5
            sigmaxE = pow(((layerWidth-i-0.5)/layerWidth),m)*sigmaMax;
            kxE = 1+(kmax-1)*pow(((layerWidth-i-0.5)/layerWidth),m);
            //parameters for c5,6H: no x shift
            sigmaxH = pow(((layerWidth-i)/layerWidth),m)*sigmaMax;
            kxH = 1+(kmax-1)*pow(((layerWidth-i)/layerWidth),m);
            //c5,6E
            Fields->c5E[i] = 2*kxE+sigmaxE*dt;
            Fields->c6E[i] = 2*kxE-sigmaxE*dt;
            //c5,6H
            Fields->c5H[i] = 2*kxH+sigmaxH*dt;
            Fields->c6H[i] = 2*kxH-sigmaxH*dt;
        }
        //simulation area
        if(layerWidth<=i && i<(Grid->Nx)-layerWidth) {
            sigmaxE = 0;
            kxE = 1;
            sigmaxH = 0;
            kxH = 1;
            Fields->c5E[i] = 2*kxE+sigmaxE*dt;
            Fields->c6E[i] = 2*kxE-sigmaxE*dt;
            Fields->c5H[i] = 2*kxH+sigmaxH*dt;
            Fields->c6H[i] = 2*kxH-sigmaxH*dt;
        }
        //right layer
        if((Grid->Nx-layerWidth-1)<=i && i<Grid->Nx) {
            sigmaxE = pow(((i+0.5-(Grid->Ny-layerWidth)+1)/layerWidth),m)*sigmaMax;
            kxE = 1+(kmax-1)*pow(((i+0.5-(Grid->Ny-layerWidth)+1)/layerWidth),m);
            sigmaxH = pow(((i-(Grid->Ny-layerWidth)+1)/layerWidth),m)*sigmaMax;
            kxH = 1+(kmax-1)*pow(((i+0.5-(Grid->Ny-layerWidth)+1)/layerWidth),m);
            Fields->c5E[i] = 2*kxE+sigmaxE*dt;
            Fields->c6E[i] = 2*kxE-sigmaxE*dt;
            Fields->c5H[i] = 2*kxH+sigmaxH*dt;
            Fields->c6H[i] = 2*kxH-sigmaxH*dt;
        }
    }
}



void Contract4x4MatrixWithFourVec(double result[4], double matrix[4][4], double vec[4])
{
    for (int i = 0; i < 4; i++)
    {
        result[i] = standard4DScalarProduct(matrix[i], vec);
    }
}

void ScaleFourVec(double a, double vec[4])
{
    int i;
    for (i = 0; i < 4; i++)
    {
        vec[i] *= a;
    }
}

double MinkowskiProduct(const double vec1[4], const double vec2[4])
{
    int i;
    double result = 0;
    for (i = 1; i < 4; i++)
    {
        result += vec1[i] * vec2[i];
    }
    result = (vec1[0] * vec2[0]) - result;
    return result;
}


void FourVecAdd(double result[4], const double vec1[4], const double vec2[4])
{
    int i;
    for (i = 0; i < 4; i++)
    {
        result[i] = vec1[i] + vec2[i];
    }
    
}

void CalcFForNystrom(double f[4], double x[4], double u[4], const double q_m, double *E, double *B, double F[4][4])
{
    Contract4x4MatrixWithFourVec(f, F, u);
}



void VecProduct(double result[3], const double vec1[3], const double vec2[3])
{
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

double ScalarProduct(const double vec1[3], const double vec2[3])
{
    int i;
    double result = 0;
    for (i = 0; i < 3; i++)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}


bool isInsideBackwardLightcone(double xParticle[4], double xObserver[4]) {
    double dt;
    double dxsq;
    double xObserverMinusxParticle[4];
    vectorDifference(xObserver, xParticle, xObserverMinusxParticle);
    dt = xObserverMinusxParticle[0];
    
    dxsq = vectorProduct(xObserverMinusxParticle, xObserverMinusxParticle);
    
    return (dt > 0 && dt*dt > dxsq);
}


bool isInsideForwardLightcone(double xParticle[4], double xObserver[4]) {
    
    double dt;
    double dxsq;
    double xObserverMinusxParticle[4];
    
    vectorDifference(xObserver, xParticle, xObserverMinusxParticle);
    dt = xObserverMinusxParticle[0];
    
    dxsq = vectorProduct(xObserverMinusxParticle, xObserverMinusxParticle);
    
    return (dt < 0 && dt*dt > dxsq);
}


void calculateIntersectionPointforTest(double xInside[4], double xOutside[4], double uInside[4], double uOutside[4], double xObserver[4], double intersectionPoint[4], double velocityAtIntersectionPoint[4]) {
    double lambda = calculateLambdaForLinearInterpolationforTest(xInside, xOutside, xObserver);
    for(int i = 0; i < 4; i++){
        intersectionPoint[i] = xOutside[i] + lambda * (xInside[i] - xOutside[i]);
        velocityAtIntersectionPoint[i] = uOutside[i] + lambda * (uInside[i] - uOutside[i]);
    }
    
}


void calculateBetaforTest(double xOld[4], double xNew[4], double beta[3]) {
    double dt = xNew[0]-xOld[0];
    beta[0] = (xNew[1]-xOld[1])/dt;
    beta[1] = (xNew[2]-xOld[2])/dt;
    beta[2] = (xNew[3]-xOld[3])/dt;
}


double calculateLambdaForLinearInterpolationforTest(double xInside[4], double xOutside[4], double xObserver[4]) {
    
    double a, b, c;
    double lambda;
    double xInsideMinusxOutside[4];
    double xOutsideMinusxObserver[4];
    
    vectorDifference(xInside, xOutside, xInsideMinusxOutside);
    vectorDifference(xOutside, xObserver, xOutsideMinusxObserver);
    
    a = minkowskiProduct(xInsideMinusxOutside, xInsideMinusxOutside);
    b = 2.0 * minkowskiProduct(xInsideMinusxOutside, xOutsideMinusxObserver);
    c = minkowskiProduct(xOutside, xOutside) + minkowskiProduct(xObserver, xObserver) - 2.0 * minkowskiProduct(xOutside, xObserver);
    
    lambda = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    return lambda;
}



void calculateLienardWiechertParametersforTest(double xParticle[4], double xObserver[4], double u[4], double *gamma_sq, double *R_sq, double *R, double n[3], double beta[3]) {
    n[0] = xObserver[1]-xParticle[1];
    n[1] = xObserver[2]-xParticle[2];
    n[2] = xObserver[3]-xParticle[3];
    
    *R_sq = n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
    *R = sqrt(*R_sq);
    
    n[0] /= *R;
    n[1] /= *R;
    n[2] /= *R;
    
    *gamma_sq = u[0]*u[0];
    
    beta[0] = u[1]/u[0];
    beta[1] = u[2]/u[0];
    beta[2] = u[3]/u[0];
}


void calculateBetaDotforTest(double uOld[4], double uNew[4], double dt, double betaDot[3]) {
    betaDot[0] = (uNew[1]/uNew[0]-uOld[1]/uOld[0])/dt;
    betaDot[1] = (uNew[2]/uNew[0]-uOld[2]/uOld[0])/dt;
    betaDot[2] = (uNew[3]/uNew[0]-uOld[3]/uOld[0])/dt;
}



void calcuateLienardWiechertFieldsforTest(double gamma_sq, double R_sq, double R, double *n,
                                   double *beta, double *beta_dot, double charge, double *E,
                                   double *B) {
    double denominator1;
    double denominator2;
    double oneMinusBetaN;
    double oneMinuBetaNCubed;
    double nMinusBeta[3];
    double nCrossnMinusBetaCrossBetaDot[3];
    double nMinusBetaCrossBetaDot[3];
    
    
    if (R_sq==0||R==0) {
        E[0] = 0;
        E[1] = 0;
        E[2] = 0;
        B[0] = 0;
        B[1] = 0;
        B[2] = 0;
        return;
    }
    
    oneMinusBetaN = 1.0-(beta[0]*n[0]+beta[1]*n[1]+beta[2]*n[2]);
    oneMinuBetaNCubed = oneMinusBetaN*oneMinusBetaN*oneMinusBetaN;
    denominator1 = 1.0/(gamma_sq*oneMinuBetaNCubed*R_sq);
    denominator2 = 1.0/(oneMinuBetaNCubed*R);
    nMinusBeta[0] = n[0]-beta[0];
    nMinusBeta[1] = n[1]-beta[1];
    nMinusBeta[2] = n[2]-beta[2];
    
    crossProduct(nMinusBeta, beta_dot, nMinusBetaCrossBetaDot);
    crossProduct(n, nMinusBetaCrossBetaDot, nCrossnMinusBetaCrossBetaDot);
    
    E[0] = charge*(denominator1*nMinusBeta[0]+denominator2*nCrossnMinusBetaCrossBetaDot[0]);
    E[1] = charge*(denominator1*nMinusBeta[1]+denominator2*nCrossnMinusBetaCrossBetaDot[1]);
    E[2] = charge*(denominator1*nMinusBeta[2]+denominator2*nCrossnMinusBetaCrossBetaDot[2]);
    
    crossProduct(n, E, B);
}


void vectorDifference(double xObserver[4], double xParticle[4], double xObserverMinusxParticle[4]) {
    xObserverMinusxParticle[0] = xObserver[0] - xParticle[0];
    xObserverMinusxParticle[1] = xObserver[1] - xParticle[1];
    xObserverMinusxParticle[2] = xObserver[2] - xParticle[2];
    xObserverMinusxParticle[3] = xObserver[3] - xParticle[3];
}


double vectorProduct(double x[4], double y[4]){
    
    return x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
}


void updateNearField(Grid *Grid, Fields *Fields, particle *particle, double time) {
    
    if (particle->didParticleChangeBoxAfterPush == true){
        printf("updating NearField ...\n");
        
        for(int i = 0; i < 27; i++){
            if(particle->boxIndicesOfNearFieldBoxesBeforePush[i] == -1 || particle->boxIndicesOfNearFieldBoxesAfterPush[i] == -1){
                continue;
            }
            if(!boxIsInNearFieldOfParticle(Grid, particle, particle->boxIndicesOfNearFieldBoxesBeforePush[i])) {
                addLWFieldsInBoxforTest(Grid, Fields, particle, particle->boxIndicesOfNearFieldBoxesBeforePush[i], time);
            }
            
            bool wasInNFBefore = false;
            for (int j = 0; j < 27; j++){
                if ( particle->boxIndicesOfNearFieldBoxesAfterPush[i] == particle->boxIndicesOfNearFieldBoxesBeforePush[j] )
                    wasInNFBefore = true;
            }
            if (wasInNFBefore == false){
                subLWFieldsInBoxforTest(Grid, Fields, particle, particle->boxIndicesOfNearFieldBoxesAfterPush[i], time);
            }
        }
    }
    
}


///@brief adds LW fields in box specified by input parameter "boxIndex"
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@param boxIndex specified box in which fields shall be calcualted
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@remark translate box index into number of box in x,y and z direction first (ib, jb, kb). Then calculate the gridIndex of the lower left corner of the current box and start looping through all grid points in that box. Each gridpoint is the current observation point, where fields shall be calcualted. Also grid index is updated each time the current grid point changes in the loop. Then add Lw field at that point.
void addLWFieldsInBoxforTest(Grid *Grid, Fields *Fields, particle *Particle, int boxIndex, double t) {
    //printf("adding LW fields in box %d\n", boxIndex);
    
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->dx;
    double dy = Grid->dy;
    double dz = Grid->dz;
    
    int adjustmentDueToUpmlInXRight = 0;
    int adjustmentDueToUpmlInXLeft = 0;
    int adjustmentDueToUpmlInYBack = 0;
    int adjustmentDueToUpmlInYFront = 0;
    int adjustmentDueToUpmlInZTop = 0;
    int adjustmentDueToUpmlInZBottom = 0;
    
    int ib = boxIndex / (numberOfBoxesInY * numberOfBoxesInZ);
    int jb = (boxIndex - ib * (numberOfBoxesInY * numberOfBoxesInZ)) / numberOfBoxesInY;
    int kb = (boxIndex - ib * (numberOfBoxesInZ * numberOfBoxesInZ)) - jb * numberOfBoxesInZ;
    
    double xObserver[4] = {0};
    xObserver[0] = t;
    
    int lowerLeftGridIndexInBox = ib * numberOfGridPointsForBoxInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + jb * numberOfGridPointsForBoxInY * numberOfGridPointsInZ * 3 + kb * numberOfGridPointsForBoxInZ * 3;
    
    for (int id = 0 + adjustmentDueToUpmlInXLeft; id < numberOfGridPointsForBoxInX - adjustmentDueToUpmlInXRight; id++){
        for (int jd = 0 + adjustmentDueToUpmlInYFront; jd < numberOfGridPointsForBoxInY - adjustmentDueToUpmlInYBack; jd++){
            for (int kd = 0 + adjustmentDueToUpmlInZBottom; kd < numberOfGridPointsForBoxInZ - adjustmentDueToUpmlInZTop; kd++){
                
                int gridIndexInBox = lowerLeftGridIndexInBox + 3 * kd + 3 * jd * numberOfGridPointsInZ + 3 * id * numberOfGridPointsInZ * numberOfGridPointsInY;
                
                xObserver[1] = (ib * numberOfGridPointsForBoxInX + id) * dx;
                xObserver[2] = (jb * numberOfGridPointsForBoxInY + jd) * dy;
                xObserver[3] = (kb * numberOfGridPointsForBoxInZ + kd) * dz;
                
                addLWField(Grid, Particle, &Fields->magneticField[gridIndexInBox], xObserver, 3);
                addLWField(Grid, Particle, &Fields->magneticField[gridIndexInBox + 1], xObserver, 4);
                addLWField(Grid, Particle, &Fields->magneticField[gridIndexInBox + 2], xObserver, 5);
                
                addLWField(Grid, Particle, &Fields->electricField[gridIndexInBox], xObserver, 0);
                addLWField(Grid, Particle, &Fields->electricField[gridIndexInBox + 1], xObserver, 1);
                addLWField(Grid, Particle, &Fields->electricField[gridIndexInBox + 2], xObserver, 2);
            }
        }
    }
}

void setInitialConditionsToStartAtT(particle *particles, int numberOfParticles) {
    for(int p = 0; p < numberOfParticles; p++) {
        for(int i = 0; i < 3; i++) {
            particles[p].xRel[i + 1] *= -1;
            particles[p].iterationCount = 0;
            particles[p].historyLength = 0;
        }
    }
}



///@brief adds LW fields in box specified by input parameter "boxIndex"
///@param Grid pointer to an instance of a Grid struct.
///@param Particle pointer to an instance of a Particle struct
///@param boxIndex specified box in which fields shall be calcualted
///@param t time at which the fields shall be calcualted. This is also the zeroth component of xObserver
///@remark translate box index into number of box in x,y and z direction first (ib, jb, kb). Then calculate the gridIndex of the lower left corner of the current box and start looping through all grid points in that box. Each gridpoint is the current observation point, where fields shall be calcualted. Also grid index is updated each time the current grid point changes in the loop. Then add Lw field at that point.
void subLWFieldsInBoxforTest(Grid *Grid, Fields *Fields, particle *Particle, int boxIndex, double t) {
    //printf("adding LW fields in box %d\n", boxIndex);
    
    
    int numberOfGridPointsForBoxInX = Grid->numberOfGridPointsPerBoxInX;
    int numberOfGridPointsForBoxInY = Grid->numberOfGridPointsPerBoxInY;
    int numberOfGridPointsForBoxInZ = Grid->numberOfGridPointsPerBoxInZ;
    
    int numberOfGridPointsInY = Grid->Ny;
    int numberOfGridPointsInZ = Grid->Nz;
    
    int numberOfBoxesInY = Grid->numberOfBoxesInY;
    int numberOfBoxesInZ = Grid->numberOfBoxesInZ;
    
    double dx = Grid->dx;
    double dy = Grid->dy;
    double dz = Grid->dz;
    
    int adjustmentDueToUpmlInXRight = 0;
    int adjustmentDueToUpmlInXLeft = 0;
    int adjustmentDueToUpmlInYBack = 0;
    int adjustmentDueToUpmlInYFront = 0;
    int adjustmentDueToUpmlInZTop = 0;
    int adjustmentDueToUpmlInZBottom = 0;
    
    int ib = boxIndex / (numberOfBoxesInY * numberOfBoxesInZ);
    int jb = (boxIndex - ib * (numberOfBoxesInY * numberOfBoxesInZ)) / numberOfBoxesInY;
    int kb = (boxIndex - ib * (numberOfBoxesInZ * numberOfBoxesInZ)) - jb * numberOfBoxesInZ;
    
    double xObserver[4] = {0};
    xObserver[0] = t;
    
    int lowerLeftGridIndexInBox = ib * numberOfGridPointsForBoxInX * numberOfGridPointsInY * numberOfGridPointsInZ * 3 + jb * numberOfGridPointsForBoxInY * numberOfGridPointsInZ * 3 + kb * numberOfGridPointsForBoxInZ * 3;
    
    for (int id = 0 + adjustmentDueToUpmlInXLeft; id < numberOfGridPointsForBoxInX - adjustmentDueToUpmlInXRight; id++){
        for (int jd = 0 + adjustmentDueToUpmlInYFront; jd < numberOfGridPointsForBoxInY - adjustmentDueToUpmlInYBack; jd++){
            for (int kd = 0 + adjustmentDueToUpmlInZBottom; kd < numberOfGridPointsForBoxInZ - adjustmentDueToUpmlInZTop; kd++){
                
                int gridIndexInBox = lowerLeftGridIndexInBox + 3 * kd + 3 * jd * numberOfGridPointsInZ + 3 * id * numberOfGridPointsInZ * numberOfGridPointsInY;
                
                xObserver[1] = (ib * numberOfGridPointsForBoxInX + id) * dx;
                xObserver[2] = (jb * numberOfGridPointsForBoxInY + jd) * dy;
                xObserver[3] = (kb * numberOfGridPointsForBoxInZ + kd) * dz;
                
                subLWField(Grid, Particle, &Fields->magneticField[gridIndexInBox], xObserver, 3);
                subLWField(Grid, Particle, &Fields->magneticField[gridIndexInBox + 1], xObserver, 4);
                subLWField(Grid, Particle, &Fields->magneticField[gridIndexInBox + 2], xObserver, 5);
                
                subLWField(Grid, Particle, &Fields->electricField[gridIndexInBox], xObserver, 0);
                subLWField(Grid, Particle, &Fields->electricField[gridIndexInBox + 1], xObserver, 1);
                subLWField(Grid, Particle, &Fields->electricField[gridIndexInBox + 2], xObserver, 2);
                
            }
        }
    }
}


///@brief calculates gamma from spatial components of given velocity vector via gamma = sqrt(1 + u^2)
///@remark Don't use this method before spatial components of u are initialized
double getGammaFromVelocityVector(double u[4]){
    double result = 0.0;
    for (int i = 1; i < 4; i++){
        result += u[i] * u[i];
    }
    return sqrt(1 + result);
    
}



