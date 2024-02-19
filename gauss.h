
extern int gauss(double **matrix, int rows, int cols);
extern double determinant(double **matrix, int rows, int cols);
extern void backward(double **UMatrix, int rowsNumber);
extern void columnSwap(double **matrix, int rowsNumber, int col1, int col2);
extern int rank(double **matrix, int dim);
extern void transpose(double **matrix, int dim);
extern double **newMatrix(int rows, int cols);
extern void matrixByVector(double **matrix, double v[], int dim);
