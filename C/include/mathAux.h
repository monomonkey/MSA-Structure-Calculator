#ifndef MATHAUX_H
#define MATHAUX_H

typedef struct {

    size_t mats_number;
    size_t mats_row_size;
    size_t mats_col_size;

    double*** tensor;

} tensor3d;

double* doubleCallocArray(int len, char* identifier);
double** allocateMatrix(int rows, int cols);
tensor3d* allocateTensor3(size_t mats_row_size, size_t mats_col_size, size_t mats_number);
void freeMatrix(double** matrix, int rows);
void freeTensor3(tensor3d* tensor);
double sumArray(double* array, size_t N);
void set_identity(double** matrix, int n);
void matrix_copy(double** dest, double** src, int rows, int cols);
void sumMatrix(double** A, double** B, double** C, size_t m, size_t n);
void diffMatrix(double** A, double** B, double** C, size_t m, size_t n);
void diagInverse(double** A, double** C, size_t m);
void elementWiseMulMatrix(double** A, double** B, double** C, size_t m, size_t n);
void multMatrix(double** A, double** B, double** C, size_t m, size_t k, size_t n);
void setVector2zeros(double* V, size_t m);
void scaleVector(double* V, double k, size_t m);
void scaleMatrix(double** A, double k, size_t m, size_t n);
void flatMatrix(double** A, double* V, size_t m, size_t n);
void unflatMatrix(double** A, double* V, size_t m, size_t n);
void inverseMatrix(double** A, size_t m, size_t n);
// vector* solveLinearSystemGSL(matrix* A, vector* x, vector* b);
void solveLinearSystemLapack(double* A, double* x, double* b, int n);
double* getDiagonal(double** matrix, size_t N);
double** buildDiagMat(double* diag, size_t N);
void printArrayd(double* v, int len, char* name);
void printArraye(double* v, int len, char* name);
void printArrayi(int* v, int len, char* name);
void printMatrixd(double** v, int len1, int len2, char* name);
// ############################################################
double FT(double k, int power, double A, double B);
double FTexp(double k, double kappa, double A, double B);
double FT(double k, int power, double A, double B);
double FTexp(double k, double kappa, double A, double B);
double FT_TaylorSeries_at_k0(double k, int power, double A, double B);
double FTexp_TaylorSeries_at_k0(double k, double kappa, double A, double B);

#endif // MATHAUX_H