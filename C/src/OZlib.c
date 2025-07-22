#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

#include <lapacke.h>
#include <cblas.h>

#include "../include/mathAux.h"
#include "../include/sysVar.h"
#include "../include/fileManagement.h"
#include "../include/OZlib.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// Check the electroneutrality of the system
bool IsElectroneutral(SysVar sysVar) {
    
    double totalSystemCharge = 0.0;
    
    for (size_t i = 0; i < sysVar.nSpecies; i++) {
        totalSystemCharge += sysVar.z[i] * sysVar.rho[i];
    }
    
    bool isElectroneutral = fabs(totalSystemCharge) < 1e-15;
    
    printf("IS THE SYSTEM ELECTRONEUTRAL?   %s\n", isElectroneutral ? "YES" : "NOPE");
    
    return isElectroneutral;
}

void rescaleRho(double phi, SysVar* sysVar) {
    
    double phi0 = 0.0;
    double proportion;

    // Calculate phiVec and sum to get phi0
    for (size_t i = 0; i < sysVar->nSpecies; i++) {

        phi0 += (M_PI / 6.0) * sysVar->rho[i] * 1e-30 * 
                (sysVar->sigma[i] * sysVar->sigma[i] * sysVar->sigma[i]);
    }

    // Calculate proportion
    proportion = phi / phi0;

    // Calculate and store rescaled rho values
    for (size_t i = 0; i < sysVar->nSpecies; i++) {
        sysVar->rho[i] = sysVar->rho[i] * proportion;
        sysVar->rhoMol[i] = sysVar->rho[i] / (Av * 1e3);

    }

    return;
}

double getPhiTotal(int nSpecies, double phi[]) {

    double phiTotal = 0.0;

    for (int i=0; i<nSpecies; i++){
        phiTotal += phi[i];
    }

    return phiTotal;

}

void phi2rho(double phi[], SysVar* sysVar) {
    
    // Calculate and store rho values
    for (size_t i = 0; i < sysVar->nSpecies; i++) {
        sysVar->rho[i] = (6.0 / M_PI) * phi[i] / (sysVar->sigma[i] * sysVar->sigma[i] * sysVar->sigma[i]);
        sysVar->rho[i] *= 1e30;

        sysVar->rhoMol[i] = sysVar->rho[i] / (Av * 1e3);
    }

    return;
}

void rescaleT(double Tstar, SysVar* sysVar) {

    double sigmaRef = sysVar->sigma[sysVar->nSpecies - 1];

    sysVar->T = e0*e0 * Tstar;
    sysVar->T /= (kB * sysVar->E * 4.0*M_PI * sigmaRef * 1e-10);

    sysVar->BETA = 1.0 / (kB * sysVar->T);

    return;
}

void rescaleE(double Tstar, SysVar* sysVar) {

    double zRef = sysVar->z[sysVar->nSpecies - 1];
    double sigmaRef = sysVar->sigma[sysVar->nSpecies - 1];

    sysVar->E = Tstar * sysVar->BETA * (zRef * zRef) * (e0*e0);
    sysVar->E /= (4.0*M_PI * sigmaRef * 1e-10);

    return;
}

void getKappaD(SysVar* sysVar) {

    // CALCULATING THE DEBYE-HUCKEL SCREENING CONSTANT

    sysVar->kappa = 0.0;

    for (size_t i = 0; i < sysVar->nSpecies; i++) {
        sysVar->kappa += e0*e0 * sysVar->z[i]*sysVar->z[i] * sysVar->rho[i];
    }

    sysVar->kappa = 1e-10 * sqrt(sysVar->BETA * sysVar->kappa / sysVar->E);

    return;
}

void BuildXi(MSAVar* msaVar, SysVar sysVar) {

    double eta[sysVar.nSpecies];

    for (size_t i = 0; i < sysVar.nSpecies; i++) {
        eta[i] = (M_PI/6.0) * sysVar.rho[i];
    }

    for (size_t i = 0; i < 4; i++) {
        msaVar->xi[i] = 0.0;

        for (size_t j = 0; j < sysVar.nSpecies; j++) {
            msaVar->xi[i] = msaVar->xi[i] + eta[j] * 
                            pow(sysVar.sigma[j]*1e-10, i);
        }
    }

    return;
}

void GammaFunction(double gammaSeed, MSAVar* msaVar, SysVar sysVar) {
    
    // Allocate memory for R
    double R[sysVar.nSpecies];

    // Cast Angstroms to meters
    for (int i = 0; i < sysVar.nSpecies; i++) {
        R[i] = sysVar.sigma[i] * 1e-10;
    }

    // Gamma Functions
    double cTilde = (M_PI / 2.0) / (1.0 - msaVar->xi[3]);
    double chiNumerator = 0.0;
    double chiDenominator = 0.0;

    for (int i = 0; i < sysVar.nSpecies; i++) {
        double zita = 1.0 + gammaSeed * R[i];
        chiNumerator += sysVar.rho[i] * R[i] * sysVar.z[i] / zita;
        chiDenominator += sysVar.rho[i] * R[i]*R[i]*R[i] / zita;
    }

    chiNumerator *= -cTilde;
    chiDenominator = 1.0 + cTilde*chiDenominator;
    double chi = chiNumerator / chiDenominator;

    // printf("chi = %.10e\n", chi);

    double DGamma = 0.0;
    for (int i = 0; i < sysVar.nSpecies; i++) {
        double zita = 1.0 + gammaSeed * R[i];
        double X = (sysVar.z[i] + R[i] * R[i] * chi) / zita;
        DGamma += sysVar.rho[i] * X * X;

        // printf("X[%d] = %.10e ", i, X);        
    }

    double omega = (sysVar.BETA * e0 * e0) / (4.0 * M_PI * sysVar.E);
    msaVar->FGamma = gammaSeed * gammaSeed - M_PI * omega * DGamma;

    // printf("\nDGamma     = %.10e\n", DGamma);
    // printf("omega      = %.10e\n", omega);
    // printf("FGamma     = %.10e\n", msaVar->FGamma);


    // Gamma Functions Derivatives
    double dchiNumerator = 0.0;
    double dchiDenominator = 0.0;

    for (int i = 0; i < sysVar.nSpecies; i++) {
        double zita = 1.0 + gammaSeed * R[i];
        dchiNumerator += sysVar.rho[i] * R[i]*R[i] * sysVar.z[i] / (zita*zita);
        dchiDenominator += sysVar.rho[i] * pow(R[i], 4.0) / (zita * zita);
    }

    dchiNumerator *= cTilde;
    dchiDenominator *= -cTilde;

    double dchi = (dchiNumerator*chiDenominator - chiNumerator*dchiDenominator) / (chiDenominator*chiDenominator);

    double dDGamma = 0.0;
    for (int i = 0; i < sysVar.nSpecies; i++) {
        double zita = 1.0 + gammaSeed * R[i];
        double X = (sysVar.z[i] + R[i] * R[i] * chi) / zita;
        double dX = (R[i]*dchi*zita - sysVar.z[i] - R[i]*R[i]*chi) * R[i] / (zita*zita);
        dDGamma += sysVar.rho[i] * 2.0 * X * dX;
    }

    msaVar->dFGamma = 2.0*gammaSeed - M_PI * omega * dDGamma;

    return;
}

void GammaNewtonRaphson(double maxError, MSAVar* msaVar, SysVar sysVar) {

    double gammaSeed = 1.0;    

    msaVar->gammaError = 1.0;

    // printf("\n----------------------------\n");
    // printf("CALCULANDO EL VALOR DE GAMMA\n");
    // printf("----------------------------\n");
    
    while (msaVar->gammaError > maxError) {
        
        GammaFunction(gammaSeed, msaVar, sysVar);

        // printf("FGamma = %.10e\n", msaVar->FGamma);
        // printf("dFGamma = %.10e\n", msaVar->dFGamma);
        
        msaVar->gamma = gammaSeed - (msaVar->FGamma / msaVar->dFGamma);
        
        msaVar->gammaError = fabs((gammaSeed - msaVar->gamma) / msaVar->gamma);

        gammaSeed = msaVar->gamma;
    }

    return;
}

double* BuildLinspace(double xFrom, double xTo, int nSteps) {
    
    if (nSteps <= 0) {
        printf("Error: Number of steps must be positive.\n");
        exit(EXIT_FAILURE);
    }

    double dx = (xTo - xFrom) / (nSteps-1);
    double* vector = (double*)malloc(nSteps * sizeof(double));

    if (vector == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < nSteps; i++) {
        vector[i] = xFrom + dx*(i);
    }

    return vector;
}

void BuildCijMatrices(int i, int j, double** cFConst, MSAVar msaVar, SysVar sysVar) {
    
    // **********CAST ANGSTROMS TO METERS**********
    double R[sysVar.nSpecies];
    double eta[sysVar.nSpecies];
    double zita[sysVar.nSpecies];

    #pragma omp parallel for
    for (int k = 0; k<sysVar.nSpecies; k++) {
        R[k] = sysVar.sigma[k] * 1e-10;
        eta[k] = (M_PI / 6.0) * sysVar.rho[k];
        zita[k] = 1.0 + msaVar.gamma * R[k];
    }
    
    // **********BASIC FUNCTIONS**********
    double Rij = (R[i] + R[j]) / 2.0;
    double lamij = fabs(R[j] - R[i]) / 2.0;
    double omega = (sysVar.BETA * e0*e0) / (4.0 * M_PI * sysVar.E);
    double cTilde = (M_PI / 2.0) / (1.0 - msaVar.xi[3]);

    double chiNumerator = 0.0;
    double chiDenominator = 0.0;

    #pragma omp parallel for reduction(+:chiNumerator, chiDenominator)
    for (int k = 0; k<sysVar.nSpecies; k++) {
        chiNumerator += sysVar.rho[k] * R[k] * sysVar.z[k] / zita[k];
        chiDenominator += sysVar.rho[k] * R[k]*R[k]*R[k] / zita[k];
    }

    chiNumerator *= -cTilde;
    chiDenominator = 1.0 + cTilde*chiDenominator;
    double chi = chiNumerator / chiDenominator;

    double X[sysVar.nSpecies];
    double N[sysVar.nSpecies];
    
    #pragma omp parallel for
    for (int k = 0; k<sysVar.nSpecies; k++) {
        X[k] = (sysVar.z[k] + R[k]*R[k] * chi) / zita[k];
        N[k] = (R[k] * chi - msaVar.gamma * sysVar.z[k]) / zita[k];
    }


    // **********Cij FUNCTIONS**********
    double beta0;
    if (R[j] > R[i]) {
        beta0 = 2.0*omega * (sysVar.z[i]*N[j] - X[i]*R[i]*chi + R[i]*R[i]*R[i] * chi*chi/3.0);
    } else {
        beta0 = 2.0*omega * (sysVar.z[j]*N[i] - X[j]*R[j]*chi + R[j]*R[j]*R[j] * chi*chi/3.0);
    }
    
    double alpha0 = omega * lamij*lamij * ((X[i] + X[j])*chi - Rij*Rij*chi*chi + N[i]*N[j]);
    
    double alpha1 = omega * ((chi*chi/3.0) * (pow(R[i], 3.0) + pow(R[j], 3.0)) -
                             (X[i] - X[j]) * (N[i] - N[j]) -
                             (X[i]*X[i] + X[j]*X[j]) * msaVar.gamma -
                             2.0 * Rij * N[i]*N[j]);
    
    double alpha2 = omega * ((X[i] + X[j])*chi + N[i]*N[j] -
                             (chi*chi/2.0) * (R[i]*R[i] + R[j]*R[j]));
    
    double alpha3 = omega * chi * chi / 3.0;
    
    double ni = 1.0 / (1.0 - msaVar.xi[3]);
    
    double a = 3.0 * msaVar.xi[2] * ni*ni;
    
    double b = 3.0 * ni*ni * (msaVar.xi[1] + (3.0 * msaVar.xi[2] * msaVar.xi[2] * ni));
    
    double c = 3.0 * ni*ni * (msaVar.xi[0] + (6.0 * msaVar.xi[1] * msaVar.xi[2] * ni) +
                              (9.0 * msaVar.xi[2]*msaVar.xi[2]*msaVar.xi[2] * ni*ni));

    double q[sysVar.nSpecies];
    double v[sysVar.nSpecies];
    #pragma omp parallel for
    for (int k = 0; k<sysVar.nSpecies; k++) {
        q[k] = ni + (a*R[k]) + (b*R[k]*R[k]) + (c * R[k]*R[k]*R[k] / 3.0);
        v[k] = (-a/2.0) - (b*R[k]) - (c * R[k]*R[k]/2.0);
    }

    // **********OTHER CONSTANTS**********
    double Uij = sysVar.z[i] * sysVar.z[j] * e0 * e0 / (4.0 * M_PI * sysVar.E);
    double pij = -(lamij*lamij/2.0) * (a + 2.0*b*Rij + c*Rij*Rij);
    double qij = (q[i] + q[j]) / 2.0;
    double vij = (v[i] + v[j]) / 2.0;

    double w = 0.0;
    #pragma omp parallel for reduction(+:w)
    for (int k = 0; k<sysVar.nSpecies; k++) {
        w += eta[k] * q[k];
    }
    w *= (1.0/2.0);

    // **********CREATION OF MATRICES**********
    for (size_t row = 0; row < 3; row++) {
        for (size_t col = 0; col < 4; col++) {
            cFConst[row][col] = 0.0;
        }
    }

    // **********FILLING THE MATRICES**********
    // filling for 0 <= s <= lamij
    if (R[j] > R[i]) {
        cFConst[0][1] = beta0 - q[i];
    } else {
        cFConst[0][1] = beta0 - q[j];
    }

    // filling for lamij < s <= Rij
    cFConst[1][0] = alpha0 - pij;
    cFConst[1][1] = alpha1 - qij;
    cFConst[1][2] = alpha2 - vij;
    cFConst[1][3] = alpha3 - w;

    // filling for s > Rij
    cFConst[2][0] = -sysVar.BETA * Uij;

    // **********CAST METERS TO ANGSTROMS**********
    for (int row = 0; row < 3; row++) {
        cFConst[row][0] *= 1e10;
        cFConst[row][2] *= 1e-10;
        cFConst[row][3] *= 1e-30;
    }

    return;
}

void writeArray(char* filename, double* array, MainVar mainVar, SysVar sysVar) {

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not open hr matrix output file :'( . \n");
        return;
        // return EXIT_FAILURE;
    }
    
    // Write the matrix to the file
    for (int i = 0; i < mainVar.totalNodes; i++) {

        fprintf(file, "%.10e ", mainVar.x[i]);
        fprintf(file, "%.10e ", array[i]);
        fprintf(file, "\n");

    }

    fclose(file);

    return;
}


// double getdPdRho(MainVar mainVar, SysVar sysVar) {

//     size_t N1 = mainVar.totalNodes;
//     size_t N2 = mainVar.nSpecies;

//     double molSum = sumArray(sysVar.rhoMol, N2);

//     double** rhoMat = allocateMatrix(N2, N2);
//     set_identity(rhoMat, N2);

//     for (size_t i=0; i<N2; i++) {
//         rhoMat[i][i] = sysVar.rhoMol[i];
//     }

//     double suma = 0.0;
//     double* crVec = doubleCallocArray(N1, "for crVec in getdPdRho");
//     double** cr = allocateMatrix(N2, N2);

//     for (size_t o=0; o<N1; o++) {

//         // Copy crMat[o] into cr
//         matrix_copy(cr, mainVar.crMat->tensor[o], N2, N2);

//         multMatrix(rhoMat, cr, cr, N2, N2, N2);
//         multMatrix(cr, rhoMat, cr, N2, N2, N2);

//         suma = 0.0;
//         for (size_t i=0; i<N2; i++) {
//             for (size_t j=0; j<N2; j++) {
//                 suma += cr[i][j];
//             }
//         }

//         crVec[o] = suma * mainVar.x[o]*mainVar.x[o];
//     }

//     double xFrom = mainVar.x[0];
//     double xTo = mainVar.x[mainVar.totalNodes-1];

//     double dPdRho = IntegralSimps(xFrom, xTo, mainVar.x, crVec, mainVar.totalNodes);
//     dPdRho = 1.0 - 1.0/(molSum) * 4.0*M_PI * dPdRho * 1e-30 * Av * 1e3;


//     free(crVec); crVec=NULL;
//     freeMatrix(rhoMat, N2);
//     freeMatrix(cr, N2);

//     return dPdRho;
// }


void BuildCkMSA(double kVec[], int kLen, double rmax, SysVar sysVar, MSAVar* msaVar, tensor3d* Cij_out) {
    
    // **********CAST ANGSTROMS TO METERS**********
    double R[sysVar.nSpecies];

    #pragma omp parallel for
    for (int k = 0; k<sysVar.nSpecies; k++) {
        R[k] = sysVar.sigma[k];// * 1e-10;
    }
    
    // Calcular parámetros xi
    BuildXi(msaVar, sysVar);

    // for(int i=0; i<4; i++) {
    //     printf("Xi[%d] = %.5e\n", i, msaVar->xi[i]);
    // }

    // Calcular gamma usando Newton-Raphson
    double gammaError = 1e-08;
    GammaNewtonRaphson(gammaError, msaVar, sysVar);
    // printf("gamma = %.10e\n\n", msaVar->gamma);

    for (int i = 0; i < sysVar.nSpecies; i++) {
        for (int j = 0; j < sysVar.nSpecies; j++) {

            // **********BASIC FUNCTIONS**********
            double Rij = (R[i] + R[j]) / 2.0;
            double lamij = fabs(R[j] - R[i]) / 2.0;

            LimitsVar limitsVar; init_limitsVar(&limitsVar);
            BuildCijMatrices(i, j, limitsVar.cFConst, *msaVar, sysVar);

            // printMatrixd(limitsVar.cFConst, 3, 4, "cFConst\n");

            for (int k = 0; k < kLen; k++) {
                double kval = kVec[k];
                double val = 0.0;

                if (kval < 0.002) {
                    val += limitsVar.cFConst[1][1] * FT_TaylorSeries_at_k0(kval, 0, lamij, Rij);
                    val += limitsVar.cFConst[1][2] * FT_TaylorSeries_at_k0(kval, 1, lamij, Rij);
                    val += limitsVar.cFConst[1][3] * FT_TaylorSeries_at_k0(kval, 3, lamij, Rij);
                    val += limitsVar.cFConst[0][1] * FT_TaylorSeries_at_k0(kval, 0, 0.0, lamij);

                    if (fabs(lamij) > 1e-10)
                        val += limitsVar.cFConst[1][0] * FT_TaylorSeries_at_k0(kval, -1, lamij, Rij);

                    if (fabs(Rij) > 0.0)
                        val += limitsVar.cFConst[2][0] * FT_TaylorSeries_at_k0(kval, -1, Rij, rmax);

                } else {
                    val += limitsVar.cFConst[1][1] * FT(kval, 0, lamij, Rij);
                    val += limitsVar.cFConst[1][2] * FT(kval, 1, lamij, Rij);
                    val += limitsVar.cFConst[1][3] * FT(kval, 3, lamij, Rij);
                    val += limitsVar.cFConst[0][1] * FT(kval, 0, 0.0, lamij);

                    if (fabs(lamij) > 1e-10)
                        val += limitsVar.cFConst[1][0] * FT(kval, -1, lamij, Rij);

                    if (fabs(Rij) > 0.0)
                        val += limitsVar.cFConst[2][0] * FT(kval, -1, Rij, rmax);
                }

                Cij_out->tensor[k][i][j] = 4.0 * M_PI * val;

                // printf("Cij = %.5e\n", val);
            }

            free_limitsVar(&limitsVar);            
        }
    }

}


void BuildCrMSA(double r[], int rLen, SysVar sysVar, MSAVar* msaVar, tensor3d* Cij_out) {

    // **********CAST ANGSTROMS TO METERS**********
    double R[sysVar.nSpecies];

    #pragma omp parallel for
    for (int k = 0; k<sysVar.nSpecies; k++) {
        R[k] = sysVar.sigma[k];// * 1e-10;
    }
    
    // Calcular parámetros xi
    BuildXi(msaVar, sysVar);

    // for(int i=0; i<4; i++) {
    //     printf("Xi[%d] = %.5e\n", i, msaVar->xi[i]);
    // }

    // Calcular gamma usando Newton-Raphson
    double gammaError = 1e-08;
    GammaNewtonRaphson(gammaError, msaVar, sysVar);
    // printf("gamma = %.10e\n\n", msaVar->gamma);

    for (int i = 0; i < sysVar.nSpecies; i++) {
        for (int j = 0; j < sysVar.nSpecies; j++) {

            // **********BASIC FUNCTIONS**********
            double Rij = (R[i] + R[j]) / 2.0;
            double lamij = fabs(R[j] - R[i]) / 2.0;

            LimitsVar limitsVar; init_limitsVar(&limitsVar);
            BuildCijMatrices(i, j, limitsVar.cFConst, *msaVar, sysVar);

            // printMatrixd(limitsVar.cFConst, 3, 4, "cFConst\n");

            for (int k = 0; k < rLen; k++) {
                double rk = r[k];
                double val = 0.0;

                if (rk < lamij) {
                    val = limitsVar.cFConst[0][1];

                } else if ((rk >= lamij) && (rk < Rij)) {
                    val = limitsVar.cFConst[1][1] +
                          limitsVar.cFConst[1][2] * rk +
                          limitsVar.cFConst[1][3] * pow(rk, 3);

                    if (fabs(lamij) > 0.0)
                        val += limitsVar.cFConst[1][0] / rk;
                } else {
                
                    if ((rk >= Rij) && (fabs(Rij) > 0.0)) 
                        val = limitsVar.cFConst[2][0] / rk;
                }

                Cij_out->tensor[k][i][j] = val;

                // printf("Cij = %.5e\n", val);
            }

            free_limitsVar(&limitsVar);
        }
    }

}
