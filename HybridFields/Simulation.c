//
//  Simulation.c
//  HybridFields
//
//  Created by Julian Manke on 13.01.17.
//  Copyright © 2017 Julian Manke. All rights reserved.
//

#include "Simulation.h"
#include "Calculations.h"
#include "Particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Grid.h"
#include "Fields.h"

//==============================================================================================================
//==============================================================================================================
//Liener-Wichert: For each time step, particle gets pushed by Nyström and send fields (Liener Wiechert), which will be calculated
//on xy-plane (for fixed z) and updated each time step
//==============================================================================================================
//==============================================================================================================


/////@brief actual simulation. Simulates: Nyström-Push, Electric-and Magnetic Fields, induced by particle movement,...etc.
//void LienerWiechert() {
//    double tN=10; //Standard Simulation Time = 10
//    //double tStart=0;
//    //double h=(tN-tStart)/number_of_Iterations_N_global;
//    
//    int numberOfParticles = 1;
//    
//    //creataes particle struct
//    //particle particle;
//    particle particles[numberOfParticles];
//    
//    //creates Grid struct
//    Grid Grid;
//    
//    //creates Fields struct
//    Fields Fields;
//    
//    
//    //initialize Grid with resolution (dx,dy,dz) and Number of Grid points (Nx,Ny,Nz).
//    //Attention not to make Nx,Ny,Nz too big in order to stay within int long range. Grid has Grid Points from (including) 0 up to Nidi in each direction
//    initGrid(&Grid,0.1, 0.1, 0.1, 16, 16, 16, 16, 16, 16);
//    double dt = 0.5*Grid.dx;
//    int numberOfIterations = tN/dt;
//    int nz = Grid.Nz;
//    int ny = Grid.Ny;
//    
//    //initialize particle with: mass m
//    initParticle(&particles[0], 1, numberOfIterations);
//    
//    //initialize Fields
//    initFields(&Fields,&Grid);
//    
//    //Initial values for DGL --> Outsource this to init function!
//    
//    //t0
//    particles[0].xRel[0] = 0;
//    //x0
//    particles[0].xRel[1] = 17;
//    //y0
//    particles[0].xRel[2] = 11;
//    //z0
//    particles[0].xRel[3] = 14.8;
//    //x0'
//    particles[0].uRel[1] = 0.458;
//    //y0'
//    particles[0].uRel[2] = 0;
//    //z0'
//    particles[0].uRel[3] = 0;
//    
//    particles[0].xRelHistory[0][0] = particles[0].xRel[0];
//    particles[0].xRelHistory[0][1] = particles[0].xRel[1];
//    particles[0].xRelHistory[0][2] = particles[0].xRel[2];
//    particles[0].xRelHistory[0][3] = particles[0].xRel[3];
//    
//    particles[0].uRelHistory[0][0] = sqrt(1+pow(particles[0].xRel[1],2)+pow(particles[0].xRel[2],2)+pow(particles[0].xRel[3],2));
//    particles[0].uRelHistory[0][1] = particles[0].uRel[1];
//    particles[0].uRelHistory[0][2] = particles[0].uRel[2];
//    particles[0].uRelHistory[0][3] = particles[0].uRel[3];
//    
//    //initialize Parameters for "calculateLienerWiechertParameters"
//    double n[3]={0};
//    double R;
//    double beta[3]={0};
//    double gamma_sq;
//    double E[3]={0};
//    double B[3]={0};
//    double betaDot[3]={0};
//    double intersectionPoint[4]={0};
//    double velocityAtIntersectionPoint[4]={0};
//    
//    //Open file to write movement of particle within Nyström Method
//    FILE *fid=fopen("Particle.txt","w");
//    
//    //time iteration ("number_of_Iterations_N_global" time steps)
//    for(int p=0;p<numberOfIterations;p++) {
//        //printf("p=%d\n",p);
//        
//        printf("Calculation of time step %d...\n", p);
//        //Nyström calculates movement for particle CHECKED: CALCULATES CORRECTLY
//        Nystrom(particles, &Fields, &Grid, p,fid,dt, 1, 0, dt, tN);
//        
//        //GridIteration Count keeping track of Grid Points: for each new Grid Point being considered (000), (100),...etc., the gridIterationGrid moves up by one and gets resetted after each Grid Point for ONE time iteration, i.e. moves along matrix columns for each row
//        //int gridIterationCount=0;
//        
//        //for Loops for all Grid Points, with xobserve[0]=x, xobserve[1]=y, xobserve[2]=z, loops like: (000), (100),...(Nx*dy00),(010),...(0 Ny*dy 0),(001),...,(00 Nz*dz)
//        //for(int p=0;p<=Grid.Nz;p++) {
//        for(int i=0;i<=Grid.Nx;i++) {
//            for(int j=0;j<=Grid.Ny;j++) {
//                
//                //redefine xobserve
//                Grid.xobserve[0]=particles[0].xRel[0]; //Test if this is properly working
//                //printf("xobserve[0]=time=%f\n",particle.actualPhase[0]);
//                Grid.xobserve[1]=i*(Grid.dx);//x-component of xobserve
//                Grid.xobserve[2]=j*(Grid.dy); //y-component of xobserve
//                Grid.xobserve[3]=14.8; //Take more general value for this! (Calculated Value)
//                //printf("xobserve[x][y][z]=%f %f %f\n",Grid.xobserve[0], Grid.xobserve[1], Grid.xobserve[2]);
//                
//                //Calculate Retarded time.
//                int retardedTime = calculationOfRetardedTime(&particles[0],Grid.xobserve);
//                
//                //Only continue calculating Parameters and fields, if retardedTime !=0, because if it is 0, i.e. that retardedTime does not exist
//                if(retardedTime!=0) {
//                    //Calculate interpolated intersectiom point with light cone
//                    calculateIntersectionPoint(retardedTime, &particles[0], Grid.xobserve, intersectionPoint, velocityAtIntersectionPoint);
//                    //Calculation of electric-and magnetic field for each time-step
//                    calculateLienerWiechertParameters(intersectionPoint, velocityAtIntersectionPoint, Grid.xobserve, n, &R, beta, &gamma_sq);
//                    //printf("xobserve[0][1][2]=%f %f %f\n",Grid.xobserve[0],Grid.xobserve[1],Grid.xobserve[2]);
//                    calculateBetaDot(&particles[0], betaDot, retardedTime, dt);
//                    calcuateLienardWiechertFields(gamma_sq, R, n, beta, betaDot, particles[0].charge, E, B); //CHECK
//                    //printf("electricField= %f %f %f\n",Fields.electricField[0],Fields.electricField[1],Fields.electricField[2]);
//                    
//                    //Save calculated Liener-Wiechert Fields in array
////                    for(int j=0;j<3;j++) {
////                        Fields.electricField[(gridIterationCount*3)+j]=E[j];
////                        //printf("%f %f %f\n",E[0],E[1],E[2]);
////                    }
//                    Fields.electricField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (int) (Grid.xobserve[2]/Grid.dz) + 0]=E[0];
//                    //printf("GridIndex=%d\n",3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (int) (Grid.xobserve[2]/Grid.dz) + 0);
//                    Fields.electricField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (int) (Grid.xobserve[2]/Grid.dz) + 1]=E[1];
//                    Fields.electricField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (int) (Grid.xobserve[2]/Grid.dz) + 2]=E[2];
//                    //printf("electric field x=%f\n",Fields.electricField[0]);
//                    
//                    //Calculate absolute value (double) of (3D) Liener-Wiechert Fields and save it in array "absoluteElectricField"
//                    //Fields.absoluteElectricField[gridIterationCount]= calculateAbsolutValueSquaredOf3DVector(E);
//                    
//                } else if(retardedTime == 0) {
//                    //fill electric field with zeros, if there is no electric field (no retarded time)
////                    for(int j=0;j<3;j++) {
////                        Fields.electricField[(gridIterationCount*3)+j]=0;
////                    }
//                    Fields.electricField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (int) (Grid.xobserve[2]/Grid.dz) + 0]=0;
//                    Fields.electricField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (int) (Grid.xobserve[2]/Grid.dz) + 1]=0;
//                    Fields.electricField[3 * (nz) * (ny) * i + 3 * (nz) * j + 3 * (int) (Grid.xobserve[2]/Grid.dz) + 2]=0;
//                    
//                    
//                    //Fields.absoluteElectricField[gridIterationCount] = 0;
//                }
//                //printf("gridIteration=%d\n",gridIterationCount);
//                //gridIterationCount=gridIterationCount+1; //Loops correctly
//            }//end of for loop for GridCoordinate-x
//        }
//        //}
//        //Before going to the next time iteration stepm write calculated fields (for each grid point) in a file (for fixed simulation time)
//        writeElectricFieldToFile(&Grid, &particles[0], &Fields, p);
//        //Iteration Count will keep track of the actual simulation (time) step, which is needed within "calculationOfRetarededTime" func
//        particles[0].iterationCount=particles[0].iterationCount+1;
//    }
//    fclose(fid);
//    
//    
//    //Test Prints
//    //printf("Electric Field component E[0][1]=%f\n",Fields.electricField[200][1]);
//    
//    //Bash Scripts:
//    //Copies FILE fid = particle movement to working path and starts python Matlib Plot to show results
//    //system("cd ~/Documents/UNI/LMU_München/Masterarbeit/Numerics/groundwork/Nystoem_Particle_in_B_Field;. BashScript.sh");
//    //Changes directory to DEBUG (XCode), opens python script for creating images and
//    system("cd ~/Documents/UNI/LMU_München/Masterarbeit/Numerics/HybridFields;. BashScript_MakeImagesAndMovieForLW.sh");
//    
//    
//    //==============================================================================================================
//    //FREE MEMORY
//    
//    
//    //deallocate Memory for electric field
//    freeMemoryforArray(Fields.electricField);
//    freeMemoryforArray(Fields.absoluteElectricField);
//    //deallocate Memory for particle y-position
//    freeMemoryforArray(particles[0].yposition);
//    //deallocate Memory for particle x-position
//    freeMemoryforArray(particles[0].xposition);
//    //deallocate Memory for particle z-position
//    freeMemoryforArray(particles[0].zposition);
//    //deallocate Memory for velocity
//    freeMemoryforArray(particles[0].xvelocity);
//    freeMemoryforArray(particles[0].yvelocity);
//    freeMemoryforArray(particles[0].zvelocity);
//    freeMemoryforArray(particles[0].velocityTime);
//    //deallocate Memory for particle time
//    freeMemoryforArray(particles[0].time);
//    //deallocate Memoty for xobserve Array
//    freeMemoryforArray(Grid.xobserve);
//    
//    
//}


//==============================================================================================================
//==============================================================================================================
//Analytic Solution: Calculates the LW fields 100% analytically
//==============================================================================================================
//==============================================================================================================

///@brief In the last calculation step, calculate all the fields analytically
void hybridFieldsWithAnalyticPrecalculation() {
    
    // =======================================
    // REQUIRED INPUT
    // =======================================
    
    //Simulation Time of actual simulation
    double tN = 15.0;
    // time of precalculation (if used)
    double precalculationTime = 0.0; // precalc of 7 works pretty well for circcles (B_z = 1)
    
    //number of particles for simulation
    int numberOfParticles = 1;
    
    char filename[32] = "";
    
    Grid Grid;
    Fields Fields;
    particle particles[numberOfParticles];
    
    // dx, numberOfBoxes, numberOfGridPointsPerBox
    initGrid(&Grid, 0.2, 0.2, 0.2, 8, 8, 8, 16, 16, 16);
    initFields(&Fields,&Grid);
    
    double dt = 0.5 * Grid.dx;
    
    int numberOfIterations = tN / dt;
    int numberOfPrecalculationSteps = precalculationTime / dt;
    
    for(int i = 0; i < numberOfParticles; i++) {
        initParticle(&particles[i], 1, numberOfIterations);
    }
    
    precalculateCoefficientsForUPML(&Grid, &Fields);
    
    // preinit all particles, on order not to get compile errors
    for(int p = 0; p < numberOfParticles; p++) {
        initializeParticle(&particles[p], 0, 0, 0, 0, 0, 0, 0);
    }
    
    // initial conditions for all particles
//    initializeParticle(&particles[1], 0, 17.0, 12.0, 14.8, -0.945, 0, 0);
    initializeParticle(&particles[0], 0, 10.5, 15.8, 14.8, 0.7, 0, 0);
    //initializeParticle(&particles[0], 0, 11.8, 11.8, 14.8, 2.5, 0, 0);
    
    // Simulation info for python script
    FILE *info = fopen("simulationInfo.txt","w");
    fprintf(info, "%f\n%f\n%f\n%i\n%i\n%i\n%i", tN, precalculationTime, Grid.dx, Grid.numberOfBoxesInX, Grid.numberOfGridPointsPerBoxInX, numberOfParticles, Grid.Nx);
    fclose(info);
    
    writeDetailedSimulationInfoToFile(particles, &Grid, numberOfParticles, precalculationTime, tN);
    
    FILE *rect = fopen("rectangleInfo.txt", "w");
    
    // actual (absolute) simulation time used for different methods
    double time = 0;
    
    // actual simulation time index
    int p=0;
    
    // define external fields
    double Bext[3] = {0};
    double Eext[3] = {0};
    
    Bext[2] = 1;
    
    // =======================================
    // ANALYTIC FIELD-PRECALCULATION
    // =======================================
    
    if(precalculationTime > 0) {
        reallocateMemoryForParticleHistories(particles, numberOfParticles, dt, tN, precalculationTime);
        nystromBackwards(particles, &Grid, &Fields, numberOfParticles, dt, precalculationTime, Bext, Eext);
        resetInitialConditions(particles, numberOfParticles, &time, Bext);
        precalculateFieldsForGivenPrecalculationTime(particles, &Grid, &Fields, numberOfParticles, numberOfPrecalculationSteps, rect, dt, tN, &p, &time, precalculationTime, Bext, Eext);
    }
    //==============================================================================================================
    // MAIN ROUTINE
    //==============================================================================================================
    for(; p < numberOfIterations + numberOfPrecalculationSteps; p++) {
        
        printf("Calculation of time step %i of %i...\n", p-numberOfPrecalculationSteps, numberOfIterations);
        writeHistoryOfParticlesToFile(particles, filename, p, numberOfParticles);
        
        pushEField(&Grid, &Fields, particles, numberOfParticles, time, dt);
        pushHField(&Grid, &Fields, particles, numberOfParticles, time + dt / 2.0, dt);
        
        for(int h=0; h < numberOfParticles; h++) {
            particles[h].oldBoxIndexBeforePush = calcCurrentBoxIndexOfParticle(&particles[h], &Grid);
            calcBoxIndizesOfNextNeighbourBoxes(&Grid, &particles[h], particles[h].boxIndicesOfNearFieldBoxesBeforePush);
            
            // write casted particle position (for rect) into file ============= --> ADJUST FOR MORE THAN TWO PARTICLES
            int numberBoxesX = (int) (particles[h].xRel[1]/(Grid.numberOfGridPointsPerBoxInX * Grid.dx));
            int numberBoxesY = (int) (particles[h].xRel[2]/(Grid.numberOfGridPointsPerBoxInY * Grid.dy));
            double cornerX = numberBoxesX * (Grid.numberOfGridPointsPerBoxInX * Grid.dx);
            double cornerY = numberBoxesY * (Grid.numberOfGridPointsPerBoxInY * Grid.dy);
            fprintf(rect, "%f %f\n", cornerX, cornerY);
            // =================================================================
//            double tau = dt / particles[h].uRel[0];
//            Nystrom(particles, &Fields, &Grid, p, tau, numberOfParticles, h, dt, tN, Bext, Eext,true);
            
//             FOR TESTING: USE BPORIS PUSHER
            updateVelocityWithBorisPusher(particles, &Grid, &Fields, numberOfParticles, h, Eext, Bext, dt, p);
            updateLocation(particles, &Grid, dt, h, p);
            
            particles[h].newBoxIndexAfterPush = calcCurrentBoxIndexOfParticle(&particles[h], &Grid);
            calcBoxIndizesOfNextNeighbourBoxes(&Grid, &particles[h], particles[h].boxIndicesOfNearFieldBoxesAfterPush);
            
            if(particles[h].oldBoxIndexBeforePush != particles[h].newBoxIndexAfterPush) {
                particles[h].didParticleChangeBoxAfterPush = true;
                updateNearField(&Grid, &Fields, &particles[h], time);
            } else {
                particles[h].didParticleChangeBoxAfterPush = false;
            }
        }
        
        pushHField(&Grid, &Fields, particles, numberOfParticles, time + dt / 2.0, dt);
        pushEField(&Grid, &Fields, particles, numberOfParticles, time, dt);
        
        writeElectricFieldToFile(&Grid, &particles[0], &Fields, p);
        
        int planeForPlotting = particles[0].xRel[3] / Grid.dz;
        writeFieldComponentsForFourierAnalysisToFile(&Grid, &Fields, filename, p, planeForPlotting, true, false);
        
        time += dt;
//        printf("p = %d\n", p);
//        printf("time = %f\n", time);
//        printf("xrel[0] = %f\n", particles[0].xRel[0]);
//        printf("xRelHistory[%d] = %f\n", p, particles[0].xRelHistory[p][0]);
//        printf("test\n");
    }
    fclose(rect);
    system("cd ~/Documents/Masterarbeit/Numerics/HybridFields;. BashScript_MakeImagesAndMovie.sh");
    
    //==============================================================================================================
    //FREE MEMORY. MAKE FUNCTION OUT OF THE FOLLOWING CODE
    for (int i = 0; i < Grid.numberOfBoxesInX * Grid.numberOfBoxesInY * Grid.numberOfBoxesInZ; i++) {
        free(Fields.Hz_im1[i]);
        Fields.Hz_im1[i] = NULL;
        free(Fields.Hy_im1[i]);
        Fields.Hy_im1[i] = NULL;
        free(Fields.Hx_jm1[i]);
        Fields.Hx_jm1[i] = NULL;
        free(Fields.Hz_jm1[i]);
        Fields.Hz_jm1[i] = NULL;
        free(Fields.Hx_km1[i]);
        Fields.Hx_km1[i] = NULL;
        free(Fields.Hy_km1[i]);
        Fields.Hy_km1[i] = NULL;
        free(Fields.Ey_ip1[i]);
        Fields.Ey_ip1[i] = NULL;
        free(Fields.Ez_ip1[i]);
        Fields.Ez_ip1[i] = NULL;
        free(Fields.Ez_jp1[i]);
        Fields.Ez_jp1[i] = NULL;
        free(Fields.Ex_jp1[i]);
        Fields.Ex_jp1[i] = NULL;
        free(Fields.Ex_kp1[i]);
        Fields.Ex_kp1[i] = NULL;
        free(Fields.Ey_kp1[i]);
        Fields.Ey_kp1[i] = NULL;
    }
    
    free(Fields.Hz_im1);
    Fields.Hz_im1 = NULL;
    free(Fields.Hy_im1);
    Fields.Hy_im1 = NULL;
    free(Fields.Hx_jm1);
    Fields.Hx_jm1 = NULL;
    free(Fields.Hz_jm1);
    Fields.Hz_jm1 = NULL;
    free(Fields.Hx_km1);
    Fields.Hx_km1 = NULL;
    free(Fields.Hy_km1);
    Fields.Hy_km1 = NULL;
    free(Fields.Ey_ip1);
    Fields.Ey_ip1 = NULL;
    free(Fields.Ez_ip1);
    Fields.Ez_ip1 = NULL;
    free(Fields.Ez_jp1);
    Fields.Ez_jp1 = NULL;
    free(Fields.Ex_jp1);
    Fields.Ex_jp1 = NULL;
    free(Fields.Ex_kp1);
    Fields.Ex_kp1 = NULL;
    free(Fields.Ey_kp1);
    Fields.Ey_kp1 = NULL;
    
    freeMemoryforArray(Fields.electricField);
    freeMemoryforArray(Fields.dField);
    freeMemoryforArray(Fields.magneticField);
    freeMemoryforArray(Fields.bField);
    freeMemoryforArray(Fields.absoluteElectricField);
    freeMemoryforArray(Fields.c1E);
    freeMemoryforArray(Fields.c2E);
    freeMemoryforArray(Fields.c3E);
    freeMemoryforArray(Fields.c4E);
    freeMemoryforArray(Fields.c5E);
    freeMemoryforArray(Fields.c6E);
    freeMemoryforArray(Fields.c1H);
    freeMemoryforArray(Fields.c2H);
    freeMemoryforArray(Fields.c3H);
    freeMemoryforArray(Fields.c4H);
    freeMemoryforArray(Fields.c5H);
    freeMemoryforArray(Fields.c6H);
    
    for(int i=0; i<numberOfParticles; i++) {
        freeMemoryforArray(particles[i].yposition);
        freeMemoryforArray(particles[i].xposition);
        freeMemoryforArray(particles[i].zposition);
        freeMemoryforArray(particles[i].xvelocity);
        freeMemoryforArray(particles[i].yvelocity);
        freeMemoryforArray(particles[i].zvelocity);
        freeMemoryforArray(particles[i].velocityTime);
        freeMemoryforArray(particles[i].time);
    }
    
    for(int j = 0;j < numberOfParticles; j++) {
        for(int i = 0; i< (precalculationTime / dt) + (tN / dt); i++) {
            free(particles[j].xRelHistory[i]);
            free(particles[j].uRelHistory[i]);
        }
        free(particles[j].xRelHistory);
        free(particles[j].uRelHistory);
    }
    //==============================================================================================================
}





//==============================================================================================================
//==============================================================================================================
//Maxwell Pusher: Tests the Maxwell Pusher for a given Pulse
//==============================================================================================================
//==============================================================================================================


//==============================================================================================================
//==============================================================================================================
//Maxwell Pusher: Tests the Maxwell Pusher for a given Pulse
//==============================================================================================================
//==============================================================================================================

///@brief ATTENTION: When choosing this method, "writingFieldsToFile()" has to be adjusted, in order to plot the right plane
void maxwellPusher() {
    
    //Simulation Time
    double tN = 10; //Standard Simulation Time = 10 --> User Input?
    
    //creates Grid struct
    Grid Grid;
    
    //creates Fields struct
    Fields Fields;
    
    //initialize Grid with resolution (dx,dy,dz) and Number of Grid points (Nx,Ny,Nz).
    //Attention not to make Nx,Ny,Nz too big in order to stay within int long range. Grid has Grid Points from (including) 0 up to Nidi in each direction
    initGrid(&Grid,0.2, 0.2, 0.2, 8, 8, 8, 16, 16, 16); //0.2, 0.2, 0.2, 8, 8, 8, 20, 20, 20 //0.125, 0.125, 0.125, 8, 8, 8, 32, 32, 32
    
    double dt = 0.5 * Grid.dx;
    int numberOfIterations = tN/dt;
    
    //creataes particle struct
    particle particle;
    
    //initialize particle with: mass m
    initParticle(&particle, 1, numberOfIterations);
    
    //initialize Fields
    initFields(&Fields,&Grid);
    
    //initialize sample pulse on grid
    initSamplePulseOnGridGausssian(&Grid, &Fields);
    //initSamplePulseOnGridSine(&Grid, &Fields);
    
    
    //time iteration ("number_of_Iterations_N_global" time steps)
    for(int p=0;p<numberOfIterations;p++) {
        printf("Calculating time step %d\n...",p);
        
        //Calculate Maxwell Pusher for the whole grid in each time iteration
        //push E
        pushEFieldOnGrid(&Fields, &Grid, dt);
        //push H
        pushHFieldOnGrid(&Fields, &Grid, dt);
        //push H
        pushHFieldOnGrid(&Fields, &Grid, dt);
        //push E
        pushEFieldOnGrid(&Fields, &Grid, dt);
        writeElectricFieldToFile(&Grid, &particle, &Fields, p);
    }
    
    
    //==============================================================================================================
    //FREE MEMORY
    
    
    //deallocate Memory for electric field
    freeMemoryforArray(Fields.electricField);
    //deallocate Memory for magnetic field
    freeMemoryforArray(Fields.magneticField);
    freeMemoryforArray(Fields.absoluteElectricField);
    //deallocate Memory for particle y-position
    freeMemoryforArray(particle.yposition);
    //deallocate Memory for particle x-position
    freeMemoryforArray(particle.xposition);
    //deallocate Memory for particle z-position
    freeMemoryforArray(particle.zposition);
    //deallocate Memory for velocity
    freeMemoryforArray(particle.xvelocity);
    freeMemoryforArray(particle.yvelocity);
    freeMemoryforArray(particle.zvelocity);
    freeMemoryforArray(particle.velocityTime);
    //deallocate Memory for particle time
    freeMemoryforArray(particle.time);
    //deallocate Memory for xobserve Array
    freeMemoryforArray(Grid.xobserve);
    
}



//==============================================================================================================
//==============================================================================================================
//==============================================================================================================
//==============================================================================================================
//Liener Wiechert Propagation with Near Field: Simulates Wave Propagation for far field 
//==============================================================================================================
//==============================================================================================================
//==============================================================================================================
//==============================================================================================================



///@brief This method simulates LW interactions within the near field area and pushed the fields via MW within the far field.
void LienerWiechertPusher() {
    
    //=============================================================================================================
    
    //Simulation Time
    double tN = 30.0;
    //number of particles for simulation
    int numberOfParticles = 2;
    
    Grid Grid;
    Fields Fields;
    particle particles[numberOfParticles];
    
    initGrid(&Grid, 0.1, 0.1, 0.1, 16, 16, 16, 16, 16, 16); //0.2, 0.2, 0.2, 8, 8, 8, 20, 20, 20 //0.125, 0.125,
    initFields(&Fields,&Grid);
    
    double dt = 0.5 * Grid.dx;
    int numberOfIterations = tN / dt;
    
    printf("number of Iterations = %i\n", numberOfIterations);
    
    for(int i = 0;i<numberOfParticles; i++) {
        initParticle(&particles[i], 1, numberOfIterations);
    }
    
    precalculateCoefficientsForUPML(&Grid, &Fields);
    
    initializeParticle(&particles[1], 0, 17.0, 12.0, 14.8, -0.958, 0, 0);
    initializeParticle(&particles[0], 0, 3.0, 11.0, 14.8, 0.958, 0, 0);
    
    //==============================================================================================================
    // actual implementation
    
    // file for particle position for Particle1
    FILE *fid1;
    fid1 = fopen("Particle1.txt", "w");
    
    // file for particle position for Particle0
    FILE *fid0;
    fid0 = fopen("Particle0.txt","w");
    
    FILE *a[2];
    a[0] = fid0;
    a[1] = fid1;
    
    // Simulation info
    FILE *info = fopen("simulationInfo.txt","w");
    fprintf(info, "%f\n%f\n%i\n%i\n%i", tN, Grid.dx, Grid.numberOfBoxesInX, Grid.numberOfGridPointsPerBoxInX, numberOfParticles);
    fclose(info);
    
    FILE *rect = fopen("rectangleInfo.txt", "w");
    
    // actual (absolute) simulation time
    double time = 0;
    
    double Bext[3] = {0};
    double Eext[3] = {0};
    
    for(int p=0; p<numberOfIterations; p++) {
        
        printf("Calculation of time step %i of %i\n", p, numberOfIterations);
        
        pushEField(&Grid, &Fields, particles, numberOfParticles, time, dt);
        pushHField(&Grid, &Fields, particles, numberOfParticles, time + dt / 2.0, dt);
        
        for(int h=0; h<numberOfParticles; h++) {
            particles[h].oldBoxIndexBeforePush = calcCurrentBoxIndexOfParticle(&particles[h], &Grid);
            calcBoxIndizesOfNextNeighbourBoxes(&Grid, &particles[h], particles[h].boxIndicesOfNearFieldBoxesBeforePush);
            
            // write casted particle position (for rect) into file ============= --> ADJUST FOR MORE THAN TWO PARTICLES
            int numberBoxesX = (int) (particles[h].xRel[1]/(Grid.numberOfGridPointsPerBoxInX * Grid.dx));
            int numberBoxesY = (int) (particles[h].xRel[2]/(Grid.numberOfGridPointsPerBoxInY * Grid.dy));
            double cornerX = numberBoxesX * (Grid.numberOfGridPointsPerBoxInX * Grid.dx);
            double cornerY = numberBoxesY * (Grid.numberOfGridPointsPerBoxInY * Grid.dy);
            fprintf(rect, "%f %f\n", cornerX, cornerY);
            // =================================================================
            double tau = dt/particles[h].uRel[0];
            Nystrom(particles, &Fields, &Grid, p, tau, numberOfParticles, h, dt, tN, Bext, Eext, true);
            
            particles[h].newBoxIndexAfterPush = calcCurrentBoxIndexOfParticle(&particles[h], &Grid);
            calcBoxIndizesOfNextNeighbourBoxes(&Grid, &particles[h], particles[h].boxIndicesOfNearFieldBoxesAfterPush);
            
            if(particles[h].oldBoxIndexBeforePush != particles[h].newBoxIndexAfterPush) {
                particles[h].didParticleChangeBoxAfterPush = true;
                updateNearField(&Grid, &Fields, &particles[h], time);
            } else {
                particles[h].didParticleChangeBoxAfterPush = false;
            }
        }
        
        pushHField(&Grid, &Fields, particles, numberOfParticles, time + dt / 2.0, dt);
        pushEField(&Grid, &Fields, particles, numberOfParticles, time, dt);
        
        writeElectricFieldToFile(&Grid, &particles[0], &Fields, p);
        
        time += dt;
    }
    
    fclose(fid0);
    fclose(fid1);
    fclose(rect);
    system("cd ~/Documents/UNI/LMU_München/Masterarbeit/Numerics/HybridFields;. BashScript_MakeImagesAndMovie.sh");
    
    //==============================================================================================================
    //FREE MEMORY. MAKE FUNCTION OUT OF THE FOLLOWING CODE
    for (int i = 0; i < Grid.numberOfBoxesInX * Grid.numberOfBoxesInY * Grid.numberOfBoxesInZ; i++) {
        free(Fields.Hz_im1[i]);
        Fields.Hz_im1[i] = NULL;
        free(Fields.Hy_im1[i]);
        Fields.Hy_im1[i] = NULL;
        free(Fields.Hx_jm1[i]);
        Fields.Hx_jm1[i] = NULL;
        free(Fields.Hz_jm1[i]);
        Fields.Hz_jm1[i] = NULL;
        free(Fields.Hx_km1[i]);
        Fields.Hx_km1[i] = NULL;
        free(Fields.Hy_km1[i]);
        Fields.Hy_km1[i] = NULL;
        free(Fields.Ey_ip1[i]);
        Fields.Ey_ip1[i] = NULL;
        free(Fields.Ez_ip1[i]);
        Fields.Ez_ip1[i] = NULL;
        free(Fields.Ez_jp1[i]);
        Fields.Ez_jp1[i] = NULL;
        free(Fields.Ex_jp1[i]);
        Fields.Ex_jp1[i] = NULL;
        free(Fields.Ex_kp1[i]);
        Fields.Ex_kp1[i] = NULL;
        free(Fields.Ey_kp1[i]);
        Fields.Ey_kp1[i] = NULL;
    }
    
    free(Fields.Hz_im1);
    Fields.Hz_im1 = NULL;
    free(Fields.Hy_im1);
    Fields.Hy_im1 = NULL;
    free(Fields.Hx_jm1);
    Fields.Hx_jm1 = NULL;
    free(Fields.Hz_jm1);
    Fields.Hz_jm1 = NULL;
    free(Fields.Hx_km1);
    Fields.Hx_km1 = NULL;
    free(Fields.Hy_km1);
    Fields.Hy_km1 = NULL;
    free(Fields.Ey_ip1);
    Fields.Ey_ip1 = NULL;
    free(Fields.Ez_ip1);
    Fields.Ez_ip1 = NULL;
    free(Fields.Ez_jp1);
    Fields.Ez_jp1 = NULL;
    free(Fields.Ex_jp1);
    Fields.Ex_jp1 = NULL;
    free(Fields.Ex_kp1);
    Fields.Ex_kp1 = NULL;
    free(Fields.Ey_kp1);
    Fields.Ey_kp1 = NULL;
    
    freeMemoryforArray(Fields.electricField);
    freeMemoryforArray(Fields.dField);
    freeMemoryforArray(Fields.magneticField);
    freeMemoryforArray(Fields.bField);
    freeMemoryforArray(Fields.absoluteElectricField);
    freeMemoryforArray(Fields.c1E);
    freeMemoryforArray(Fields.c2E);
    freeMemoryforArray(Fields.c3E);
    freeMemoryforArray(Fields.c4E);
    freeMemoryforArray(Fields.c5E);
    freeMemoryforArray(Fields.c6E);
    freeMemoryforArray(Fields.c1H);
    freeMemoryforArray(Fields.c2H);
    freeMemoryforArray(Fields.c3H);
    freeMemoryforArray(Fields.c4H);
    freeMemoryforArray(Fields.c5H);
    freeMemoryforArray(Fields.c6H);
    
    for(int i=0; i<numberOfParticles; i++) {
        freeMemoryforArray(particles[i].yposition);
        freeMemoryforArray(particles[i].xposition);
        freeMemoryforArray(particles[i].zposition);
        freeMemoryforArray(particles[i].xvelocity);
        freeMemoryforArray(particles[i].yvelocity);
        freeMemoryforArray(particles[i].zvelocity);
        freeMemoryforArray(particles[i].velocityTime);
        freeMemoryforArray(particles[i].time);
    }
    
    for(int j = 0;j < numberOfParticles; j++) {
        for(int i = 0; i< numberOfIterations; i++) {
            free(particles[j].xRelHistory[i]);
            free(particles[j].uRelHistory[i]);
        }
        free(particles[j].xRelHistory);
        free(particles[j].uRelHistory);
    }
    //==============================================================================================================
}




