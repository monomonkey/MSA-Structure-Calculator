#ifndef OZLIB_H
#define OZLIB_H

#include "sysVar.h"

bool IsElectroneutral(SysVar sysVar);
void rescaleRho(double phi, SysVar* sysVar);
double getPhiTotal(int nSpecies, double phi[]);
void phi2rho(double phi[], SysVar* sysVar);
void rescaleT(double Tstar, SysVar* sysVar);
void rescaleE(double Tstar, SysVar* sysVar);
void getKappaD(SysVar* sysVar);
void BuildXi(MSAVar* msaVar, SysVar sysVar);
void GammaFunction(double gammaSeed, MSAVar* msaVar, SysVar sysVar);
void GammaNewtonRaphson(double maxError, MSAVar* msaVar, SysVar sysVar);
double* BuildLinspace(double xFrom, double xTo, int nSteps);
void BuildCijMatrices(int i, int j, double** cFConst, MSAVar msaVar, SysVar sysVar);

void writeArray(char* filename, double* array, MainVar mainVar, SysVar sysVar);
void writeWr_KbT(char* filename, MainVar mainVar, SysVar sysVar);
void writeResidualMatrix(char* filename, MainVar mainVar, SysVar sysVar);
void writeJacobianVector(char* filename, MainVar mainVar, SysVar sysVar);
double getResidualIntegral(int j, int k, MainVar mainVar, MSAVar msaVar, SysVar sysVar);
void fillTensor(size_t centralSpecies, tensor3d* T, double** M);

double getdPdRho(MainVar mainVar, SysVar sysVar);

void BuildCkMSA(double kVec[], int kLen, double rmax, SysVar sysVar, MSAVar* msaVar, tensor3d* Cij_out);
void BuildCrMSA(double r[], int rLen, SysVar sysVar, MSAVar* msaVar, tensor3d* Cij_out);



#endif // OZLIB_H