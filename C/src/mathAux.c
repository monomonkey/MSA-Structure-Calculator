#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <lapacke.h>
#include <cblas.h>

#include "../include/mathAux.h"

double* doubleCallocArray(int len, char* identifier) {

    double* vector = (double*)calloc(len, sizeof(double));

    if (vector == NULL) {
        fprintf(stderr, "Memory allocation failed for %s.\n", identifier);
        exit(1);
    }

    return vector;
}

// Helper function to allocate memory for a 2D array
double** allocateMatrix(int rows, int cols) {

    // Allocate memory for 'rows' number of pointers (for each row)
    double** matrix = (double**) calloc(rows, sizeof(double*));

    if (matrix == NULL) {
        printf("Memory allocation failed for matrix rows\n");
        // return 1;
    }

    // Allocate memory for each row
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*) calloc(cols, sizeof(double));
        if (matrix[i] == NULL) {
            printf("Memory allocation failed for row %d\n", i);
            // return 1;
        }
    }

    return matrix;
}

// Function to allocate memory for a 3D array
tensor3d* allocateTensor3(size_t mats_row_size, size_t mats_col_size, size_t mats_number) {

    tensor3d* tensor;
    // Allocate memory for the tensor3d struct
    tensor = (tensor3d*)malloc(sizeof(tensor3d));

    // Allocate memory for "mats_number" (pointers to 2D arrays)
    tensor->tensor = (double***)malloc(mats_number * sizeof(double**));

    tensor->mats_number = mats_number;
    tensor->mats_row_size = mats_row_size;
    tensor->mats_col_size = mats_col_size;

    if (tensor->tensor == NULL) {
        perror("Memory allocation failed for tensor3 mats_number");
        return tensor;
    }

    // Allocate memory for "mats_row_size" (pointers to rows of 2D arrays)
    for (int i = 0; i < mats_number; i++) {
        tensor->tensor[i] = (double**)malloc(mats_row_size * sizeof(double*));
        if (tensor->tensor[i] == NULL) {
            perror("Memory allocation failed for tensor3 cols");
            return tensor;
        }
    
        // Allocate memory for "mats_col_size" (actual values in the 3D array)
        for (int j = 0; j < mats_row_size; j++) {
            tensor->tensor[i][j] = (double*)calloc(mats_col_size, sizeof(double));
            if (tensor->tensor[i][j] == NULL) {
                perror("Memory allocation failed for tensor3 mats_col_size");
                return tensor;
            }
        }
    }

    return tensor;
}

// Helper function to free memory for a 2D array
void freeMatrix(double** matrix, int rows) {

    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }

    free(matrix); matrix = NULL;

    return;
}

// Function to free the memory of a 3D array
void freeTensor3(tensor3d* tensor) {
    
    if (tensor->tensor == NULL) return;

    // Free the memory of mats_col_size
    for (int i = 0; i < tensor->mats_number; i++) {
        for (int j = 0; j < tensor->mats_row_size; j++) {
            free(tensor->tensor[i][j]);
        }

        // Free the memory of mats_row_size
        free(tensor->tensor[i]);
    }

    // Free the memory of mats_number
    free(tensor->tensor);

    // Free the memory of tensor
    free(tensor); tensor = NULL;

    return;
}

double sumArray(double* array, size_t N) {

    double suma = 0.0;

    for (size_t i=0; i<N; i++){
        suma += array[i];
    }

    return suma;
}

// Function to set a double** matrix to an identity matrix
void set_identity(double** matrix, int n) {

    // Initialize all elements to 0
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0.0;
        }
    }

    // Set diagonal elements to 1
    for (int i = 0; i < n; i++) {
        matrix[i][i] = 1.0;  // Set the i-th diagonal element
    }

    return;
}

// Function to copy matrix A to matrix B
void matrix_copy(double** dest, double** src, int rows, int cols) {

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            dest[i][j] = src[i][j];  // Copy element by element
        }
    }

    return;
}

/* sumMatrix performs C = A+B where A and B are mxn matrixes */
void sumMatrix(double** A, double** B, double** C, size_t m, size_t n) {

    size_t i, j;

    for (i=0; i < m; i++) {
        for (j=0; j < n; j++) {

            C[i][j] = A[i][j] + B[i][j];

        }
    }

    return;
}

/* diffMatrix performs C = A-B where A and B are mxn matrixes */
void diffMatrix(double** A, double** B, double** C, size_t m, size_t n) {

    size_t i, j;

    for (i=0; i < m; i++) {
        for (j=0; j < n; j++) {

            C[i][j] = A[i][j] - B[i][j];

        }
    }

    return;
}

/* diagInverse performs C = 1/diag(A) * ID where A is an mxn diagonal matrix */
void diagInverse(double** A, double** C, size_t m) {

    size_t i;

    matrix_copy(C, A, m, m);

    for (i=0; i < m; i++) {

        C[i][i] = 1.0 / A[i][i];

    }

    return;
}

/* elementWiseMulMatrix performs C = A*B, but element wise, where A and B are mxn matrixes */
void elementWiseMulMatrix(double** A, double** B, double** C, size_t m, size_t n) {

    size_t i, j;

    for (i=0; i < m; i++) {
        for (j=0; j < n; j++) {

            C[i][j] = A[i][j] * B[i][j];

        }
    }

    return;
}

/* multMatrix performs C = A*B where A is a mxk matrix and B is a kxn matrix */
void multMatrix(double** A, double** B, double** C, size_t m, size_t k, size_t n) {

    double* flatA = (double*)calloc(m*k, sizeof(double));
    double* flatB = (double*)calloc(k*n, sizeof(double));
    double* flatC = (double*)calloc(m*n, sizeof(double));

    if (flatA == NULL || flatB == NULL || flatC == NULL) {
        fprintf(stderr, "Memory allocation failed in mulMatrix.\n");
        exit(1);
    }

    flatMatrix(A, flatA, m, k);
    flatMatrix(B, flatB, k, n);
    flatMatrix(C, flatC, m, n);

   // Perform the matrix multiplication C = 1.0 * A * B + 0.0 * C
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k,
                1.0, flatA, k,  // A is m x k
                flatB, n,       // B is k x n
                0.0, flatC, n); // C is m x n

    unflatMatrix(C, flatC, m, n);

    free(flatA);
    free(flatB);
    free(flatC);

    return;
}

/* setVector2zeros sets all entries of aa m-sizee V to zero */
void setVector2zeros(double* V, size_t m) {

    size_t i;

    for (i=0; i < m; i++) {
            V[i] = 0.0;
    }

    return;
}

/* scaleVector performs V = kV where V is a m-size vector k is a scalar */
void scaleVector(double* V, double k, size_t m) {

    size_t i;

    for (i=0; i < m; i++) {
            V[i] = k * V[i];
    }

    return;
}

/* scaleMatrix performs A = kA where A is a mxn matrix and k is a scalar */
void scaleMatrix(double** A, double k, size_t m, size_t n) {

    size_t i, j;

    for (i=0; i < m; i++) {
        for (j=0; j < n; j++) {

            A[i][j] = k * A[i][j];

        }
    }

    return;
}

/* flatMatrix cast a mxn matrix into a m*n vector */
void flatMatrix(double** A, double* V, size_t m, size_t n) {

    size_t i, j;

    for (i=0; i < m; i++) {
        for (j=0; j < n; j++) {

            V[i*n + j] = A[i][j];

        }
    }

    return;
}

/* unflatMatrix cast a m*n vector into a mxn matrix */
void unflatMatrix(double** A, double* V, size_t m, size_t n) {

    size_t i, j;

    for (i=0; i < m; i++) {
        for (j=0; j < n; j++) {

            A[i][j] = V[i*n + j];

        }
    }

    return;
}

void inverseMatrix(double** A, size_t m, size_t n) {

    double* flatA = (double*)calloc(m*n, sizeof(double));

    flatMatrix(A, flatA, m, n);

    int* ipiv = (int*)malloc(n * sizeof(int));  // Pivot indices
    int info;  // Output info from LAPACK routines

    // Perform LU decomposition of A (dgetrf)
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, flatA, n, ipiv);
    
    if (info != 0) {
        printf("Error: LU decomposition failed, matrix may be singular\n");
        free(flatA);
        free(ipiv);

        return;
    }

    // Compute the inverse based on LU decomposition (dgetri)
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, flatA, n, ipiv);
    if (info != 0) {
        printf("Error: Matrix inversion failed\n");
    } else {
        unflatMatrix(A, flatA, m, n);
    }

    // Clean up
    free(flatA);
    free(ipiv);

    return;
}

/* This function solves the linear equation Ax = b */
// vector* solveLinearSystemGSL(matrix* A, vector* x, vector* b) {
    
//     int n = A->matrix->size1;
//     int s;  // sign of the permutation during LU decomposition

//     // Create a copy of A because gsl_linalg_LU_decomp modifies the matrix
//     // gsl_matrix* A_copy = gsl_matrix_alloc(n, n);
//     // gsl_matrix_memcpy(A_copy, A->matrix);

//     // Create a permutation vector for LU decomposition
//     gsl_permutation* p = gsl_permutation_alloc(n);

//     // Perform LU decomposition
//     gsl_linalg_LU_decomp(A->matrix, p, &s);

//     // Solve the system A x = b
//     gsl_linalg_LU_solve(A->matrix, p, b->vector, x->vector);

//     // Free memory
//     gsl_permutation_free(p);
//     // gsl_matrix_free(A_copy);

//     return x;
// }

/* This function solves the linear equation Ax = b using LAPACK 
   n is the number of elements in x and b*/
void solveLinearSystemLapack(double* A, double* x, double* b, int n) {
    
    int *ipiv = (int*)malloc(n * sizeof(int));  // Pivot indices
    int nrhs = 1;  // Number of right-hand sides (i.e., b is a vector)
    int info;

    // Copy b to x because LAPACK overwrites the solution in b.
    memcpy(x, b, n * sizeof(double));

    // Call LAPACK's dgesv function to solve the system
    // dgesv solves Ax = B using LU decomposition, with A overwritten by its factors, and B overwritten by the solution x.
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, x, nrhs);

    // Check for success
    if (info > 0) {
        printf("Singular matrix, solution not found.\n");
    }

    // Clean up
    free(ipiv);
}

// Function to extract the diagonal of a matrix
double* getDiagonal(double** matrix, size_t N) {

    double* diag = doubleCallocArray(N, "for diag in getDiagonal");
    
    for (size_t i=0; i<N; i++) {
        diag[i] = matrix[i][i];
    }

    return diag;
}

// Function to create a diagonal matrix from an array
double** buildDiagMat(double* diag, size_t N) {
    
    double** diagMat = allocateMatrix(N, N);

    for (size_t i=0; i<N; i++) {
        diagMat[i][i] = diag[i];  // Set diagonal element
    }
    
    return diagMat;
}

void printArrayd(double* v, int len, char* name) {

    printf("%s", name);
    printf("[");
    
    for (int i=0; i<len; i++){
        printf("%.5lf, ", v[i]);
    }
    
    printf("\b\b]\n");

    return;
}

void printArraye(double* v, int len, char* name) {

    printf("%s", name);
    printf("[");
    
    for (int i=0; i<len; i++){
        printf("%.10e, ", v[i]);
    }
    
    printf("\b\b]\n");

    return;
}

void printArrayi(int* v, int len, char* name) {

    printf("%s", name);
    printf("[");
    
    for (int i=0; i<len; i++){
        printf("%d, ", v[i]);
    }
    
    printf("\b\b]\n");

    return;
}

void printMatrixd(double** v, int len1, int len2, char* name) {

    printf("%s", name);
    printf("[");
    
    for (int i=0; i<len1; i++){

        (i == 0) ? printf("[") : printf(" [");

        for (int j=0; j<len2; j++){
            printf("%.5e, ", v[i][j]);
        }

        (i == (len1-1)) ? printf("\b\b]]\n") : printf("\b\b],\n");
    }

    return;
}


// ################################################################################
// FÓRMULAS ANALÍTICAS Y MÉTODOS NUMÉRICOS
// ################################################################################

double FT(double k, int power, double A, double B) {
    
    double result = 0.0;

    if (power == -1) {
        result = (-cos(k * B) + cos(k * A)) / (k * k);
    }
    else if (power == 0) {
        result = ((sin(k * B) - k * B * cos(k * B)) -
                  (sin(k * A) - k * A * cos(k * A))) / pow(k, 3);
    }
    else if (power == 1) {
        result = (((2.0 - pow(k * B, 2)) * cos(k * B) + 2.0 * k * B * sin(k * B)) -
                  ((2.0 - pow(k * A, 2)) * cos(k * A) + 2.0 * k * A * sin(k * A))) / pow(k, 4);
    }
    else if (power == 3) {
        double kB = k * B;
        double kA = k * A;
        result = (
            (4.0 * k * B * (pow(kB, 2) - 6.0) * sin(kB) -
             (pow(kB, 4) - 12.0 * pow(kB, 2) + 24.0) * cos(kB)) -
            (4.0 * k * A * (pow(kA, 2) - 6.0) * sin(kA) -
             (pow(kA, 4) - 12.0 * pow(kA, 2) + 24.0) * cos(kA))
        ) / pow(k, 6);
    }

    return result;
}


double FTexp(double k, double kappa, double A, double B) {
    
    double denom = k * (kappa * kappa + k * k);

    double termB = exp(-kappa * B) * (kappa * sin(k * B) + k * cos(k * B));
    double termA = exp(-kappa * A) * (kappa * sin(k * A) + k * cos(k * A));

    double result = (-termB + termA) / denom;

    return result;
}

double FT_TaylorSeries_at_k0(double k, int power, double A, double B) {
    
    double result = 0.0;

    if (power == -1) {
        result = (pow(B, 2) - pow(A, 2)) / 2.0
               - pow(k, 2) * (pow(B, 4) - pow(A, 4)) / 24.0
               + pow(k, 4) * (pow(B, 6) - pow(A, 6)) / 720.0;
    }
    else if (power == 0) {
        result = (pow(B, 3) - pow(A, 3)) / 3.0
               - pow(k, 2) * (pow(B, 5) - pow(A, 5)) / 30.0
               + pow(k, 4) * (pow(B, 7) - pow(A, 7)) / 840.0;
    }
    else if (power == 1) {
        result = (pow(B, 4) - pow(A, 4)) / 4.0
               - pow(k, 2) * (pow(B, 6) - pow(A, 6)) / 36.0
               + pow(k, 4) * (pow(B, 8) - pow(A, 8)) / 960.0;
    }
    else if (power == 3) {
        result = (pow(B, 6) - pow(A, 6)) / 6.0
               - pow(k, 2) * (pow(B, 8) - pow(A, 8)) / 48.0
               + pow(k, 4) * (pow(B, 10) - pow(A, 10)) / 1200.0;
    }

    return result;
}


double FTexp_TaylorSeries_at_k0(double k, double kappa, double A, double B) {
    // Compute powers and exponentials for B
    double KB = kappa * B;
    double KB2 = KB * KB;
    double KB3 = KB2 * KB;
    double KB4 = KB3 * KB;
    double KB5 = KB4 * KB;
    double expKB = exp(-KB);

    // Compute powers and exponentials for A
    double KA = kappa * A;
    double KA2 = KA * KA;
    double KA3 = KA2 * KA;
    double KA4 = KA3 * KA;
    double KA5 = KA4 * KA;
    double expKA = exp(-KA);

    // First order term
    double AB = -expKB * (KB + 1.0) / (kappa * kappa);
    double AA = -expKA * (KA + 1.0) / (kappa * kappa);

    // Second order term (k^2)
    double BB = k * k * expKB * (KB3 + 3.0 * KB2 + 6.0 * KB + 6.0) / (6.0 * pow(kappa, 4));
    double BA = k * k * expKA * (KA3 + 3.0 * KA2 + 6.0 * KA + 6.0) / (6.0 * pow(kappa, 4));

    // Fourth order term (k^4)
    double k4 = k * k * k * k;
    double CB = k4 * expKB * (KB5 + 5.0 * KB4 + 20.0 * KB3 + 60.0 * KB2 + 120.0 * KB + 120.0) / (120.0 * pow(kappa, 6));
    double CA = k4 * expKA * (KA5 + 5.0 * KA4 + 20.0 * KA3 + 60.0 * KA2 + 120.0 * KA + 120.0) / (120.0 * pow(kappa, 6));

    // Final result
    double result = (AB - AA) + (BB - BA) - (CB - CA);

    return result;
}
