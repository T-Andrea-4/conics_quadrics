//
//  diagonalization.c
//  conics_quadrics
//
//  Created by Andrea Turchetto on 04/02/24.
//

#include "diagonalization.h"
#include <math.h>
#define MAX_SWEEPS 20

void rotation(double **matrix, int dim, double **eigenVectors, int r, int c);
int sgn(double d);
void jacobi(double **matrix, double **diagonalMatrix, double **eigenVectors, int dim);

int sgn(double d)
{
    if(d >= 0)
    {
        return 1;
    }

    return -1;
}

void jacobi(double **matrix, double **diagonalMatrix, double **eigenVectors, int dim) /*dim is the number of rows and columns of the matrix. This function modify the diagonalMatrix pointer (prevously allocated as a dynamic matrix)*/
{
    
    int i, j, h, diagonal; /*diagonal is the flag that stores whether  the matrix is diagonal or not */
    
    /*initialization of the eigenvectors matrix that initially is the identity matrix*/
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            if(i == j)
            {
                eigenVectors[i][j] = 1;
            }
            else
            {
                eigenVectors[i][j] = 0;
            }
        }
    }
    

    
    /*copying matrix in diagonalMatrix*/
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            diagonalMatrix[i][j] = matrix[i][j];
        }
    }

    diagonal = 0;
    for(h = 0; h < 50 && !diagonal; h++){
        
        diagonal = 1;
        for(i = 0; i < dim; i++)
        {
            for(j = i + 1; j < dim; j++)
            {
                if(diagonalMatrix[i][j] != 0.0)
                {
                    diagonal = 0;
                    rotation(diagonalMatrix, dim, eigenVectors, i, j);

                }
            }
        }
    }
}


void rotation(double **matrix, int dim, double **eigenVectors, int r, int c)
{
    double t, theta, cos, sin, tau, temp;
    int i;

    /*finding the angle needed to set the off diagonal element m[r][c] to 0*/
    theta = (matrix[r][r] - matrix[c][c]) / (2 * matrix[r][c]);
    
    
    t = sgn(theta) / (fabs(theta) + sqrt(pow(theta,2) + 1));
    

    cos = 1 / sqrt(pow(t,2) + 1);
    sin = t * cos;

    /*modifying the off diagonal elements and the diagonal elements with the values of sin and cos found before*/
    matrix[r][r] = matrix[r][r] + t * matrix[c][r];
    matrix[c][c] = matrix[c][c] - t * matrix[c][r];
    matrix[r][c] = 0.0;
    matrix[c][r] = 0.0;

    tau = sin / (1 + cos);

    for(i = 0; i < dim; i++)
    {
        temp = matrix[i][c];
        if(i != c  && i != r)
        {
            matrix[i][c] = matrix[i][c] - sin * (matrix[i][r] + (tau * matrix[i][c]));
            matrix[c][i] = matrix[i][c];
            matrix[i][r] = matrix[i][r] + sin * (temp - (tau * matrix[i][r]));
            matrix[r][i] = matrix[i][r];
        }
    }


    //modifying eigenvectors
    for(i = 0; i < dim; i++)
    {
        temp = eigenVectors[i][c];
        eigenVectors[i][c] = cos * eigenVectors[i][c] - sin * eigenVectors[i][r];
        eigenVectors[i][r] = sin * temp  + cos * eigenVectors[i][r];
    }
}
