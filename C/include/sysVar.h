#ifndef SYSVAR_H
#define SYSVAR_H

#include <stdbool.h>

#include "mathAux.h"

// Define the structure
typedef struct {
    int nSpecies;

    double T;
    double Er;
    double E;
    double BETA;
    double kappa;

    double* z;      // valence
    double* rho;    // numeric density (#/m^3)
    double* rhoMol; // molarity
    double* sigma;  // diameter    
} SysVar;

typedef struct {
    int totalNodes;
    int headSegments;
    int tailSegments;
    int nodesPerSegment;
    int* headNodesVec;
    int* tailNodesVec;
    
    double xTailMax;
    double xTailProportion;
    double* tailPointsVec;
    
    bool explicitNodes;
    bool explicitBreaks;    
} MeshVar;

typedef struct {

    double gamma;
    double gammaError;
    double FGamma;
    double dFGamma;

    double* xi;

} MSAVar;

typedef struct {

    double** integralLimits;
    double** analytLimits;
    double** cFConst;

} LimitsVar;

typedef struct {

    int maxCounter;
    double maxError;

    int totalNodes;
    int nSpecies;
    int totalEfectiveNodes;
    int* efectiveNodesStart;

    double* x;
    double* k;
    double** cr;
    double** hr;
    double** deltahr;
    double** residual;
    
    double* Jacobian;

    tensor3d* hrMat;
    tensor3d* crMat;
    tensor3d* SkMat;
    tensor3d* CkMat;

    double* psix;
    double* sigmax;
    double* qx;

} MainVar;

// Declaration of global constants
extern const double e0;
extern const double kB;
extern const double Av;
extern const double E0;

// Function prototypes
void init_sysVar(SysVar* sysVar, double* valence, double* molarity, double* diameter);
void free_sysVar(SysVar* sysVar);

void init_meshVar(MeshVar* meshVar, int* headNodesVec, int* tailNodesVec, double* tailPointsVec);
void free_meshVar(MeshVar* meshVar);

void init_msaVar(MSAVar* msaVar);
void free_msaVar(MSAVar* msaVar);

void init_limitsVar(LimitsVar* limitsVar);
void free_limitsVar(LimitsVar* limitsVar);

void init_mainVar(MainVar* mainVar, MeshVar meshVar);
void free_mainVar(MainVar* mainVar);

void BuildEfectiveNodes(MainVar* mainVar, MeshVar meshVar);

#endif // SYSVAR_H