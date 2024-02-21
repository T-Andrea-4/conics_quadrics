
extern int gauss(int rows, int cols, double matrix[rows][cols]);
extern double determinant( int rows, int cols, double matrix[rows][cols]);
extern void backward( int rowsNumber,int colsNumber, double UMatrix[rowsNumber][colsNumber]);
extern void columnSwap( int rowsNumber, int colsNumber, double matrix[rowsNumber][colsNumber], int col1, int col2);
extern int rank( int dim, double matrix[dim][dim]);
extern void transpose( int dim, double matrix[dim][dim]);
void matrixByVector( double v[], int dim, double matrix[dim][dim]);
