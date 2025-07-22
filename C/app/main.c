#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>

#include "../include/sysVar.h"
#include "../include/OZlib.h"
#include "../include/fileManagement.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

int main() {

    // SET PARAMETERS FOR PARALLEL COMPUTING -----------------------------------
    int num_procs = omp_get_num_procs();  // Get the number of available cores
    omp_set_num_threads(num_procs);       // Set OpenMP to use all cores

    //DECLARE MAIN VARIABLES ---------------------------------------------------
    SysVar sysVar;
    MSAVar msaVar;

    enum { nSpecies = 3};
    sysVar.nSpecies  = nSpecies;
    double Tstar = 9.799522328266514;
    double sigma[nSpecies] = {1, 50., 70.};
    double valence[nSpecies] = {-1.25, 8.04, 15.44};
    double phi[nSpecies] = {1.70322912e-05, 2.63762883e-01, 9.60829933e-02};
    
    double phiTotal = getPhiTotal(sysVar.nSpecies, phi);

    init_sysVar(&sysVar, valence, phi, sigma);
    init_msaVar(&msaVar);

    // CHECK SYSTEM'S ELECTRONEUTRALITY -----------------------------------------
    IsElectroneutral(sysVar);

    // BUILD MESH ---------------------------------------------------------------
    int totalNodes = 1000;
    int ktotalNodes = 100000;
    double* xVector = BuildLinspace(1e-3, 1000, totalNodes);
    double* kVector = BuildLinspace(0, 10, ktotalNodes);

    // printArraye(sysVar.rho, nSpecies, "rho ");

    // BUILD TENSORTS ---------------------------------------------------------------
    tensor3d* Cr = allocateTensor3((size_t) nSpecies, 
                                   (size_t) nSpecies, 
                                   (size_t) totalNodes);

    tensor3d* Ck = allocateTensor3((size_t) nSpecies, 
                                   (size_t) nSpecies, 
                                   (size_t) ktotalNodes);

    BuildCrMSA(xVector, totalNodes, sysVar, &msaVar, Cr);
    BuildCkMSA(kVector, ktotalNodes, xVector[totalNodes-1], sysVar, &msaVar, Ck);

    writeTensor("Cr.txt", xVector, Cr, sysVar);
    writeTensor("Ck.txt", kVector, Ck, sysVar);

    // // SET true IF WANT TO LOAD A PREVIOUS SOLUTION -----------------------------

    // init_mainVar(&mainVar, meshVar);
    // mainVar.maxCounter = maxCounter;
    // mainVar.maxError = maxError;

    // char baseFolder[100] = "OutputFiles";
    // createFolder(baseFolder);
    
    // // COMPUTING h(r) PROFILES
    // HNCMSAMultiSolver(xVector, loadSeed, hrInputfilename, mainVar, sysVar, meshVar);

    // GetMainFeatures(&mainVar, sysVar);
    // WriteMainFeatures(baseFolder, mainVar, sysVar);
    // SaveInputFilename(inputVarFile, baseFolder);

    // // Free allocated memory ########################################################
    free(xVector); xVector = NULL;
    free(kVector); kVector = NULL;
 
    free_sysVar(&sysVar);
    free_msaVar(&msaVar);

    freeTensor3(Cr);
    freeTensor3(Ck);

    return 0; 
}