#include <math.h>
#include <stdlib.h>
#include "gauss.h"

int gauss(double **matrix, int rows, int cols)
{
    int r,c, k, detSign;
    double prevPivot, firstElRow;
    int min = minimum(cols, rows);
    detSign = 1;

    for(c = 0; c < min-1 ; c++) //min-1
    {
        r = c;

        detSign *= rowSwap(matrix, r, c, rows); //moves the row with the greater pivot on top

        prevPivot = matrix[c][c]; //prevPivot is the value of the pivot in the column

        for(r = c+1; r < min; r++) //min
        {
            firstElRow = matrix[r][c];

            if(prevPivot){
                for(k = c; k < cols; k++) //k=0
                {
                    matrix[r][k] = matrix[r][k] - (matrix[c][k] * firstElRow / prevPivot);
                }
            }

        }
    }

    return detSign;

}

int minimum(int a, int b)
{
    if(b > a)
    {
        return a;
    }
    return b;
}

int rowSwap(double **matrix, int rowIndex, int colIndex, int rowsNumber)
{
    int i, rowOfMax;
    double max;
    max = matrix[rowIndex][colIndex];
    rowOfMax = rowIndex;


    //finding the greater pivot in the column
    for(i = rowIndex + 1; i < rowsNumber; i++ )
    {
        if(fabs(matrix[i][colIndex]) > fabs(max))
        {
            max = matrix[i][colIndex];
            rowOfMax = i;
        }
    }

    //swapping
    if(rowOfMax != rowIndex)
    {
        double *temp = *(matrix + rowIndex);
        *(matrix + rowIndex) = *(matrix + rowOfMax);
        *(matrix + rowOfMax) = temp;
        return -1; //returns -1 because the determinant changes sign if a row is swapped with another
    }

    return 1;


}

double determinant(double **matrix, int rows, int cols)
{
    double **Umatrix; //Umatrix is the upper matrix that results after applying gaussian elimination to matrix
    int detsign, i, j;
    double determinant;
    Umatrix = (double**) malloc(sizeof(double *) * rows);
    *Umatrix = (double *) malloc(sizeof(double ) * cols * rows);

    for(i = 0; i < rows; i++ )
    {
        *(Umatrix + i) = *(Umatrix) + cols * i;
    }

    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            Umatrix[i][j] = matrix[i][j];
        }
    }



    detsign = gauss(Umatrix, rows, cols);

    determinant = detsign;

    for (i = 0; i < rows; i++)
    {
        determinant *= Umatrix[i][i];
    }

    return determinant;

}

void backward(double **UMatrix, int rowsNumber) //backward substitution in a nx(n+1) upper triangular matrix
{
    int i,j,lastColumn;

    lastColumn = rowsNumber;

    for(i = rowsNumber - 1; i >= 0; i--)
    {
        for(j = i + 1; j < rowsNumber; j++)
        {
            UMatrix[i][lastColumn] -= -1 * UMatrix[i][j] * UMatrix[j][lastColumn];
            UMatrix[i][j] = 0;
        }

        UMatrix[i][lastColumn] /= UMatrix[i][i];
        UMatrix[i][i] = 1;
    }

}

