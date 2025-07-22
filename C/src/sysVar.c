#include <stdio.h>
#include <stdlib.h>

#include "../include/sysVar.h"
#include "../include/mathAux.h"
#include "../include/OZlib.h"

// Definition of global constants
const double e0 = 1.60217662 * 1e-19;       // Coulomb (C)
const double kB = 1.38064852 * 1e-23;       // J/K
const double Av = 6.02214129 * 1e23;        // # / mol
const double E0 = 8.8541878176 * 1e-12;     // (C^2/(J*m))


// Define functions
void init_sysVar(SysVar* sysVar, double* valence, double* phi, double* diameter) {

    sysVar->Er = 78.5;
    sysVar->T = 298.0;

    sysVar->kappa = 0.0;
    sysVar->E = E0 * sysVar->Er;
    sysVar->BETA = 1.0 / (kB * sysVar->T);

    // Allocate memory for arrays
    sysVar->z = (double*)calloc(sysVar->nSpecies, sizeof(double));
    sysVar->rho = (double*)calloc(sysVar->nSpecies, sizeof(double));
    sysVar->rhoMol = (double*)calloc(sysVar->nSpecies, sizeof(double));
    sysVar->sigma = (double*)calloc(sysVar->nSpecies, sizeof(double));

    if (sysVar->z == NULL || sysVar->rho == NULL || sysVar->rhoMol == NULL || sysVar->sigma == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1); // Exit if memory allocation fails
    }

    // Set values in the struct arrays
    for (size_t i = 0; i < sysVar->nSpecies; i++) {
        sysVar->z[i] = valence[i];
        sysVar->sigma[i] = diameter[i];
    }

    phi2rho(phi, sysVar);

    return;
    
}

void free_sysVar(SysVar* sysVar) {
    // Free allocated memory
    free(sysVar->z);
    free(sysVar->rho);
    free(sysVar->rhoMol);
    free(sysVar->sigma);

    return;
}

// ------------------------------------------------------------------------

void init_meshVar(MeshVar *meshVar, int* headNodesVec, int* tailNodesVec, double* tailPointsVec) {

    // Allocate memory for arrays based on headSegments and tailSegments
    meshVar->headNodesVec = (int*)calloc(meshVar->headSegments, sizeof(int));
    meshVar->tailNodesVec = (int*)calloc(meshVar->tailSegments, sizeof(int));
    meshVar->tailPointsVec = (double*)calloc(meshVar->tailSegments, sizeof(double));

    if (meshVar->headNodesVec == NULL || meshVar->tailNodesVec == NULL || meshVar->tailPointsVec == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1); // Exit if memory allocation fails
    }

    // Set values in the struct arrays
    for (size_t i = 0; i < meshVar->headSegments; i++) {
        meshVar->headNodesVec[i] = headNodesVec[i];
    }

    for (size_t i = 0; i < meshVar->tailSegments; i++) {
        meshVar->tailNodesVec[i] = tailNodesVec[i];
        meshVar->tailPointsVec[i] = tailPointsVec[i];
    }

    return;
}

void free_meshVar(MeshVar *meshVar) {
    // Free allocated memory
    free(meshVar->headNodesVec);
    free(meshVar->tailNodesVec);
    free(meshVar->tailPointsVec);

    return;
}

// ------------------------------------------------------------------------

void init_msaVar(MSAVar* msaVar) {
    // Allocate memory
    msaVar->xi = (double*)calloc(4, sizeof(double));

    return;
}

void free_msaVar(MSAVar* msaVar) {
    // Free allocated memory
    free(msaVar->xi);

    return;
}
// ------------------------------------------------------------------------

void init_limitsVar(LimitsVar* limitsVar) {
    // Allocate memory
    limitsVar->integralLimits = allocateMatrix(5, 2);
    limitsVar->analytLimits = allocateMatrix(5, 2);
    limitsVar->cFConst = allocateMatrix(3, 4);
    
    return;
}

void free_limitsVar(LimitsVar* limitsVar) {
    // Free allocated memory
    freeMatrix(limitsVar->integralLimits, 5);
    freeMatrix(limitsVar->analytLimits, 5);
    freeMatrix(limitsVar->cFConst, 3);

    return;
}
// ------------------------------------------------------------------------

void init_mainVar(MainVar* mainVar, MeshVar meshVar) {

    mainVar->maxCounter = 15;
    mainVar->maxError = 1e-8;

    mainVar->nSpecies = meshVar.headSegments;
    mainVar->totalNodes =  meshVar.totalNodes;
    mainVar->efectiveNodesStart = (int*)calloc(meshVar.headSegments, sizeof(int));

    BuildEfectiveNodes(mainVar, meshVar);

    mainVar->x = (double*)calloc(meshVar.totalNodes, sizeof(double));
    mainVar->k = (double*)calloc(meshVar.totalNodes, sizeof(double));

    mainVar->cr = allocateMatrix(meshVar.headSegments, meshVar.totalNodes);
    mainVar->hr = allocateMatrix(meshVar.headSegments, meshVar.totalNodes);
    mainVar->deltahr = allocateMatrix(meshVar.headSegments, meshVar.totalNodes);
    mainVar->residual = allocateMatrix(meshVar.headSegments, meshVar.totalNodes);
    
    mainVar->Jacobian = (double*)calloc(mainVar->totalEfectiveNodes * mainVar->totalEfectiveNodes, sizeof(double));

    printf("headSegments = %d\n", meshVar.headSegments);
    printf("total Nodes  = %d\n", meshVar.totalNodes);

    mainVar->hrMat = allocateTensor3((size_t) meshVar.headSegments, 
                                     (size_t) meshVar.headSegments, 
                                     (size_t) meshVar.totalNodes);
    mainVar->crMat = allocateTensor3((size_t) meshVar.headSegments, 
                                     (size_t) meshVar.headSegments, 
                                     (size_t) meshVar.totalNodes);
    mainVar->SkMat = allocateTensor3((size_t) meshVar.headSegments, 
                                     (size_t) meshVar.headSegments, 
                                     (size_t) meshVar.totalNodes);
    mainVar->CkMat = allocateTensor3((size_t) meshVar.headSegments, 
                                     (size_t) meshVar.headSegments, 
                                     (size_t) meshVar.totalNodes);

    mainVar->psix = (double*)calloc(meshVar.totalNodes, sizeof(double));
    mainVar->sigmax = (double*)calloc(meshVar.totalNodes, sizeof(double));
    mainVar->qx = (double*)calloc(meshVar.totalNodes, sizeof(double));

    return;
}

void free_mainVar(MainVar* mainVar) {

    free(mainVar->efectiveNodesStart);
    
    free(mainVar->x);
    free(mainVar->k);

    freeMatrix(mainVar->cr, mainVar->nSpecies);
    freeMatrix(mainVar->hr, mainVar->nSpecies);
    freeMatrix(mainVar->deltahr, mainVar->nSpecies);
    freeMatrix(mainVar->residual, mainVar->nSpecies);
    
    free(mainVar->Jacobian);

    freeTensor3(mainVar->hrMat);
    freeTensor3(mainVar->crMat);
    freeTensor3(mainVar->SkMat);
    freeTensor3(mainVar->CkMat);
    
    free(mainVar->psix);
    free(mainVar->sigmax);
    free(mainVar->qx);

    mainVar = NULL;

    return;
}
// ------------------------------------------------------------------------

void BuildEfectiveNodes(MainVar* mainVar, MeshVar meshVar) {
    
    mainVar->efectiveNodesStart[0] = meshVar.headNodesVec[0];
    mainVar->totalEfectiveNodes = mainVar->totalNodes - mainVar->efectiveNodesStart[0];

    for (size_t i = 1; i < meshVar.headSegments; i++) {
        mainVar->efectiveNodesStart[i] = mainVar->efectiveNodesStart[i-1] + 
                                         meshVar.headNodesVec[i];

        mainVar->totalEfectiveNodes += mainVar->totalNodes - mainVar->efectiveNodesStart[i];
    }

    return;
}