#include <math.h>
#include "gauss.h"

int rowSwap( int rowIndex, int colIndex, int rowsNumber, int colsNumber, double matrix[rowsNumber][colsNumber]);
int minimum(int a, int b);

int gauss(int rows, int cols, double matrix[rows][cols])//gauss elimination method
{
    int r,c, k, detSign;
    double prevPivot, firstElRow;
    int min = minimum(cols, rows);
    detSign = 1;

    for(c = 0; c < min-1 ; c++)
    {
        r = c;
        detSign *= rowSwap(r, c, rows, cols, matrix); //moves the row with the greater pivot on top
        prevPivot = matrix[c][c]; //prevPivot is the value of the pivot in the column c

        for(r = c + 1; r < min; r++)
        {
            firstElRow = matrix[r][c];

            if(prevPivot){
                for(k = c; k < cols; k++)
                {
                    matrix[r][k] = matrix[r][k] - (matrix[c][k] * firstElRow / prevPivot);
                }
            }
            else
            {
                for(int i = c + 1; i < cols && !prevPivot; i++)
                {
                    prevPivot = matrix[c][c];
                    firstElRow = matrix[r][c];
                }
                if(prevPivot){
                    for(k = c; k < cols; k++)
                    {
                        matrix[r][k] = matrix[r][k] - (matrix[c][k] * firstElRow / prevPivot);
                    }
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

int rowSwap( int rowIndex, int colIndex, int rowsNumber, int colsNumber, double matrix[rowsNumber][colsNumber])
{
    int i, rowOfMax;
    double max, temp;
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
        for(i = 0; i < colsNumber; i++)
        {
            temp = matrix[rowIndex][i];
            matrix[rowIndex][i] = matrix[rowOfMax][i];
            matrix[rowOfMax][i] = temp;
        }
        return -1; //returns -1 because the determinant changes sign if a row is swapped with another
    }
    return 1;
}

double determinant( int rows, int cols, double matrix[rows][cols])
{
    double Umatrix[rows][cols]; //Umatrix is the upper triangular matrix that results after applying gaussian elimination to matrix
    int detsign, i, j;
    double determinant;

    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            Umatrix[i][j] = matrix[i][j];
        }
    }

    detsign = gauss( rows, cols, Umatrix);
    determinant = detsign;

    for (i = 0; i < rows; i++)
    {
        determinant *= Umatrix[i][i];
    }

    return determinant;
}

void backward( int rowsNumber,int colsNumber, double UMatrix[rowsNumber][colsNumber]) //backward substitution in an n x (n+1) upper triangular matrix
{
    int i,j,lastColumn;
    lastColumn = rowsNumber;

    for(i = rowsNumber - 1; i >= 0; i--)
    {
        for(j = i + 1; j < rowsNumber; j++)
        {
            UMatrix[i][lastColumn] -= UMatrix[i][j] * UMatrix[j][lastColumn];
            UMatrix[i][j] = 0;
        }
        UMatrix[i][lastColumn] /= UMatrix[i][i];
        UMatrix[i][i] = 1;
    }
}

void columnSwap( int rowsNumber, int colsNumber, double matrix[rowsNumber][colsNumber], int col1, int col2)
{
    double temp;
    int i;
    
    for(i = 0; i < rowsNumber; i++)
    {
        temp = matrix[i][col1];
        matrix[i][col1] = matrix[i][col2];
        matrix[i][col2] = temp;
    }
}

int rank( int dim, double matrix[dim][dim])
{
    int found, i, j;
    double temp[dim][dim];
    
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            temp[i][j] = matrix[i][j];
        }
    }

    gauss( dim, dim, temp);

    found = 0;
    for(i = dim - 1; i >= 0 && !found; i--)
    {
        found = 0;
        for(j = dim - 1; j >= i; j--)
        {
            if(temp[i][j] != 0)
            {
                found = 1;
                i++;
            }
        }
        
    }
    
    if(found)
    {
        return i + 1;
    }
    return  0;



}

void transpose( int dim, double matrix[dim][dim])
{
    int i, j;
    double temp;
    
    for(i = 0; i < dim; i++)
    {
        for(j = i + 1; j < dim; j++)
        {
            if(j != i)
            {
                temp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = temp;
            }
        }
    }
}

void matrixByVector( double v[], int dim, double matrix[dim][dim])//computes the product rows by cols between a dim x dim matrix and a dim-long vector
{
    int i, j;
    double temp[dim];
    for(i = 0; i < dim; i++)
    {
        temp[i] = v[i];
    }

    for(i = 0; i < dim; i++)
    {
        v[i] = 0;
        for(j = 0; j < dim; j++)
        {
            v[i] += temp[j] * matrix[i][j] ;
        }
    }
}

