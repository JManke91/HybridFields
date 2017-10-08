//
//  Calculations.h
//  HybridFields
//
//  Created by Julian Manke on 27.11.16.
//  Copyright © 2016 Julian Manke. All rights reserved.
//

#ifndef Calculations_h
#define Calculations_h

#include <stdio.h>
#include <stdbool.h>
#include "Particle.h"
#include "Grid.h"
#include "Fields.h"

//extern int number_of_Iterations_N_global;
extern double deltaT;

//Definition of 2nd order DGL
//void dglfunction(int n,double input[n+1], double sol[n], double velocityTime, double E[3], double B[3], int p);

// TESTING PURPOSES
void updateVelocityWithBorisPusher(particle *Particles, Grid *Grid, Fields * Fields, int numberOfParticles, int particleIndex, double Eextern[3], double Bextern[3], double dt, int actualTimeStep);
void updateLocation(particle *Particle, Grid *Grid, double dt, int currentParticle, int currentTimeStep);

int BoxIsInNearFieldOfParticle(Grid *g, const int boxIndex, const int particleIndex);

void calcFieldsOnGridWithoutNearField(particle *particles, Grid *Grid, Fields *Fields, int numberOfParticles, double t);

void addCurrentStateToParticleHistory(particle *particle,  int index);

double E_x(const double x[4]);

double E_y(const double x[4]);

double E_z(const double x[4]);

double B_x(const double x[4]);

double B_y(const double x[4]);

double B_z(const double x[4]);

double getGammaFromVelocityVector(double u[4]);

void writeDetailedSimulationInfoToFile(particle *particles, Grid *Grid, int numberOfParticles, double precalculationTime, double calculationTime);

void writeFieldComponentsForFourierAnalysisToFile(Grid *Grid, Fields *Fields, char *filename, int index, int planeForPlotting, bool plotE, bool plotB);

void writeHistoryOfParticlesToFile(particle *particles, char *filename, int actualTimeIndex, int numberOfParticles);

void setInitialConditionsToStartAtT(particle *particles, int numberOfParticles);

void resetInitialConditions(particle *particles, int numberOfParticles, double *Bext);

void scaleVectorByConstantFactor(double vector[3], double constFactor);

void nystromBackwards(particle *particles, Grid *Grid, Fields *Fields, int numberOfParticles, double dt, double precalculationTime, double *externBField, double *externEField);

bool double_equals(double a, double b);

void setInitialConditionsToStartAtT(particle *particles, int numberOfParticles);

void reallocateMemoryForParticleHistories(particle *particles, int numberOfParticles, double dt, double simulationTime, double precalculationTime);

void precalculateFieldsForGivenPrecalculationTime(particle *particles, Grid *Grid, Fields *Fields, int numberOfParticles, int numberOfPrecalculationSteps, FILE *rect, double dt, double tN, double precalculationTime, double Bext[3], double Eext[3]);

//void updateLocation(particle *Particle, Grid *Grid, double dt, int actualTimeIndex);

double standard4DScalarProduct(double inputOne[4], double inputTwo[4]);

//void updateVelocityWithBorisPusher(particle *Particles, Grid *Grid, Fields *Fields, int numberOfParticles, int particleIndex, double dt, int actualTimeIndex);

void UpdateNNAndFields(Grid *g, particle *p, Fields *Fields, int indicesNNOld[27], int indicesNNNew[27], double t);

void updateNearField(Grid *Grid, Fields *Fields, particle *particle, double time);

void addLWFieldsInBoxforTest(Grid *Grid, Fields *Fields, particle *Particle, int boxIndex, double t);

void subLWFieldsInBoxforTest(Grid *Grid, Fields *Fields, particle *Particle, int boxIndex, double t);

void VecProduct(double result[3], const double vec1[3], const double vec2[3]);

void PushParticleVay(particle particles[], int particleIndex, int numberOfParticles,
                     Grid *g, double h, int *boxIndexBeforePush, int *boxIndexAfterPush, Fields *Fields, int currentTimeIndex, FILE *fid);

void savetyNystrom(particle *particle, Fields *Fields, Grid *Grid, int p, FILE *fid, double dt, int numberOfParticles, int actualParticle);

double ScalarProduct(const double vec1[3], const double vec2[3]);

void Contract4x4MatrixWithFourVec(double result[4], double matrix[4][4], double vec[4]);

//void CalcFForNystromLL(double f[4], const double x[4], const double u[4], const double q_m,double *E, double *B, double F[4][4]);

//void getCurrentPositionAndVelocityVector(particle *particle, double x[4], double u[4], int currentTimeIndex);

void PushParticleNystrom(particle particles[], int particleIndex, int numberOfParticles,
                        Grid *g, double ds, int *boxIndexBeforePush, int *boxIndexAfterPush, int currentTimeIndex, FILE *fid);

void UpdateFieldTensorExt(double F[4][4], double E[3], double B[3]);

//void CalcFForNystromLL(double f[4], const double x[4], const double u[4], const double q_m);

void CalcFForNystrom(double f[4], double x[4], double u[4], const double q_m, double *E, double *B, double F[4][4]);

double MinkowskiProduct(const double vec1[4], const double vec2[4]);

void FourVecAdd(double result[4], const double vec1[4], const double vec2[4]);

int GetBoxIndexForParticle(particle *p, Grid *g);

double vectorProduct(double x[4], double y[4]);

void ScaleFourVec(double a, double vec[4]);

void calcuateLienardWiechertFieldsforTest(double gamma_sq, double R_sq, double R, double *n,
                                          double *beta, double *beta_dot, double charge, double *E,
                                          double *B);

void calculateBetaDotforTest(double uOld[4], double uNew[4], double dt, double betaDot[3]);

int IsInBackwardLightCone(const double x_p[4], const double x_q[4]);

void vectorDifference(double xObserver[4], double xParticle[4], double xObserverMinusxParticle[4]);

bool isInsideForwardLightcone(double xParticle[4], double xObserver[4]);

bool isInsideBackwardLightcone(double xParticle[4], double xObserver[4]);

double calculateLambdaForLinearInterpolationforTest(double xInside[4], double xOutside[4], double xObserver[4]);

void calculateLienardWiechertParametersforTest(double xParticle[4], double xObserver[4], double u[4], double *gamma_sq, double *R_sq, double *R, double n[3], double beta[3]);

void calculateBetaforTest(double xOld[4], double xNew[4], double beta[3]);

void calculateIntersectionPointforTest(double xInside[4], double xOutside[4], double uInside[4], double uOutside[4], double xObserver[4], double intersectionPoint[4], double velocityAtIntersectionPoint[4]);

//Nyström Algorithm
void Nystrom(particle *particle, Fields *Fields, Grid *Grid, int i, double dt, int numberOfParticles, int actualParticle, double realdt, int tN, double Bext[3], double Eext[3], bool addHistoryToFile);

void analyticNystrom(particle *particle, Fields *Fields, Grid *Grid, int p, FILE *fid, double dt, int numberOfParticles, int actualParticle);

void freeMemoryforArray(double array[]);

void testPrintArray(double** array[],double arraySizeX,double arraySizeY);

int calculationOfRetardedTime(particle *particle,double xobserve[3]);

void calculateLienerWiechertParameters(double intersectionPoint[4], double velocityAtIntersectionPoint[4],double xObserve[3],double n[3],double *R,double beta[3],double *gamma_sq);

void calculateLienerWiechertParametersForCertainComponent(double intersectionPoint[4], double velocityAtIntersectionPoint[4],double xObserve[4],double n[3],double *R,double beta[3],double *gamma_sq);

void crossProduct(double a[3], double b[3], double result[3]);

void calculateLWFieldsNew(Grid *Grid, Fields *Fields, particle *particle, double dt, int component);

void calculateBetaDot(particle *particle,double betaDot[3],int retardedTimeIndex,double h);

void calcuateLienardWiechertFields(double gamma_sq, double R, double *n,
                                   double *beta, double *beta_dot, double charge, double *E,
                                   double *B);

double calculateAbsolutValueSquaredOf3DVector(double vector[3]);

void NystromNEW(particle *particle, Fields *Fields, Grid *Grid, int p, FILE *fid, double dt, int numberOfParticles, int actualParticle, double realdt, int tN);

int getCorrespondingIndexInAbsoluteElectricFieldForCoordinates(Grid *Grid, int x, int y, int z);

int getCorrespondingIndexInAbsoluteElectricFieldForCoordinatesFullGrid(Grid *Grid, int x, int y, int z);

void writeElectricFieldToFile(Grid *Grid, particle *particle,Fields *Fields, int index);

void calculateIntersectionPoint(int retardedTime,particle *particle,double xObserve[3], double intersectionPoint[4], double velocityAtIntersectionPoint[4]);

//void initSamplePulseOnGrid(Grid *Grid,Fields *Fields);

void initializeParticle(particle *particle, double t0, double x0, double y0, double z0, double xDot0, double yDot0, double zDot0);

//void initSamplePulseOnGrid2(Grid *Grid, Fields *Fields);

void initSamplePulseOnGridGausssian(Grid *Grid, Fields *Fields);

void pushEField(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, double t, double dt);

void pushHField(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, double t, double dt);


void pushEFieldOnGrid(Fields *Fields,Grid *Grid, double dt);

void pushHFieldOnGrid(Fields *Fields, Grid *Grid,double dt);

bool doesArrayContainValue(int value, int arr[], int size);

void updateFieldsAfterParticlePush(particle *particle, Fields *Fields, Grid *Grid, int nextNeighborBoxIndicesBeforePush[27], int nextNeighborBoxIndicesAfterPush[27], double dt, double time, FILE *pos);

void addLWFieldsInBox(Grid *Grid, Fields *Fields, particle *particle, int boxIndex, double time);

void setLWFieldsInBoxToZero(Grid *Grid, Fields *Fields, int boxIndex);

void subLWFieldsInBox(Grid *Grid, Fields *Fields, particle *particle, int boxIndex, double time);

void precalculateCoefficientsForUPML(Grid *Grid, Fields *Fields);

double minkowskiProduct(double x[4], double y[4]);

void pushEFieldInsideBoxes(Grid *Grid, Fields *Fields, double dt);

void pushHFieldInsideBoxes(Grid *Grid, Fields *Fields, double dt);

void setHFieldOnBorders(Grid *Grid, Fields *Fields);

void setEFieldOnBorders(Grid *Grid, Fields *Fields);

void adjustHFields(Grid *Grid, Fields *Fields, particle *particles, int numberOfParticles, const double t);

void adjustEFields(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const double t);

void adjustHyz_im1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);

void adjustHxz_jm1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);

void adjustHxy_km1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);

void adjustEyz_ip1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);

void adjustExz_jp1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);

void adjustExy_kp1(Grid *Grid, Fields *Fields, particle *particle, int numberOfParticles, const int boxIndex, const int ib, const int jb, const int kb, const double t);

void pushEFieldAtBorders(Grid *Grid, Fields *Fields, double dt);

void pushHFieldAtBorders(Grid *Grid, Fields *Fields, double dt);

int calcBoxIndexIm1(Grid *Grid, const int boxIndex);

int calcBoxIndexJm1(Grid *Grid, const int boxIndex);

int calcBoxIndexKm1(Grid *Grid, const int boxIndex);

int calcBoxIndexIp1(Grid *Grid, const int boxIndex);

int calcBoxIndexJp1(Grid *Grid, const int boxIndex);

int calcBoxIndexKp1(Grid *Grid, const int boxIndex);

bool boxIsInNearFieldOfParticle(Grid *Grid, particle *Particle, int boxIndex);

void calcBoxIndizesOfNextNeighbourBoxes(Grid *Grid, particle *particle, int boxIndizesOfNextNeighbourBoxes[27]);

int calcCurrentBoxIndexOfParticle(particle *particle, Grid *Grid);

void addLWField(Grid *Grid, particle *particle, double *destination, double xObserver[4], int component);

void subLWField(Grid *Grid, particle *particle, double *destination, double xObserver[4], int component);


#endif /* Calculations_h */
