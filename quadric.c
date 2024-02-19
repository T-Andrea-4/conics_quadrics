
#include <stdlib.h>
#include <math.h>
#include "quadric.h"
#include "gauss.h"
#include "diagonalization.h"

void quadricCanonicalForm(double **A, double **ACanonical, double **isometry, double *centre)
{
    int i,j, rkB;
    double detB, detA;
    double traslation[3];
    detB = determinant(A, 3, 3);
    detA = determinant(A, 4, 4);
    rkB = rank(A, 3);

    
    jacobi(A, ACanonical, isometry, 3);
    ACanonical[3][3] = A[3][3];

    for(i = 0; i < 4; i++) //setting the off-diagonal elements of Acanonical to 0
    {
        for(j = 0; j < 4; j++)
        {
            if(i != j)
            {
                ACanonical[i][j] = 0;
            }
        }
    }
    
    if(detB!= 0) //quadric with centre
    {
        ACanonical[3][3] = detA / (ACanonical[0][0] * ACanonical[1][1] * ACanonical[2][2]);
        
        double **centreSystem = newMatrix(3,4);

        for(i = 0; i < 3; i++)
        {
            for(j = 0; j < 4; j++)
            {
                if(j == 3)
                {
                    centreSystem[i][j] = -1 * A[i][j];
                   
                }
                else
                {
                    centreSystem[i][j] = A[i][j];
                    
                }
            }
        }
        
        gauss(centreSystem, 3, 4);
        backward(centreSystem, 3);
        centre[0] = centreSystem[0][3];
        centre[1] = centreSystem[1][3];
        centre[2] = centreSystem[2][3];
        
        free(*centreSystem);
        free(centreSystem);
        
        isometry[3][3] = 1;
        isometry[0][3] = centre[0];
        isometry[1][3] = centre[1];
        isometry[2][3] = centre[2];
        isometry[3][0] = isometry[3][1] = isometry[3][2] = 0;
    }
    else if(rkB == 2)
    {
        centre = NULL;
        /*changing the order of the eigenvalues of A so that it is simpler to find the canonical form*/
        if(ACanonical[2][2] != 0)
        {
            for(i = 0; i < 2; i++)
            {
                if(ACanonical[i][i] == 0)
                {
                    ACanonical[i][i] = ACanonical[2][2];
                    ACanonical[2][2] = 0;
                    columnSwap(isometry, 3, i, 2);
                }
            }
        }
        
        /*applying the rotation found with the jacobi method to the first grade part of A and updating the canonical matrix*/
        for(i = 0; i < 3; i++)
        {
            ACanonical[i][3] = 0;
            for(j = 0; j < 3; j++)
            {
                ACanonical[i][3] += A[j][3] * isometry[j][i];

            }

            ACanonical[3][i] = ACanonical[i][3];
        }

        /*making the translation that zeroes the terms of first grade in x and y*/
        traslation[0] = -ACanonical[0][3] / ACanonical[0][0];
        traslation[1] = -ACanonical[1][3] / ACanonical[1][1];
        ACanonical[3][3] = A[3][3] - pow(ACanonical[0][3], 2) / ACanonical[0][0] - pow(ACanonical[1][3],2) / ACanonical[1][1];
        ACanonical[0][3] = ACanonical[3][0] = 0;
        ACanonical[1][3] = ACanonical[3][1] = 0;

        if(ACanonical[2][3] != 0)
        {
            traslation[2] = -ACanonical[3][3] / (2 * ACanonical[2][3]);
            ACanonical[3][3] = 0;
        }
        matrixByVector(isometry, traslation,3);

        for(i = 0; i < 3; i++)
        {
            isometry[i][3] += traslation[i];
        }
        isometry[3][3] = 1;
    }
    else
    {
        centre = NULL;

        if(ACanonical[0][0] == 0)
        {
            for(i = 1; i < 3; i++ )
            {
                if(ACanonical[i][i] != 0)
                {
                    ACanonical[0][0] = ACanonical[i][i];
                    ACanonical[i][i] = 0;
                    columnSwap(isometry, 3, 0, i);
                }
            }
        }


        
        /*applying the rotation found with the jacobi method to the first grade part of A and storing the result in canonical matrix*/
        for(i = 0; i < 3; i++)
        {
            ACanonical[i][3] = 0;
            for(j = 0; j < 3; j++)
            {
                ACanonical[i][3] += A[j][3] * isometry[j][i];
                
            }
            
            ACanonical[3][i] = ACanonical[i][3];
        }



        /*making the traslation that zeroes the term of first grade in x*/
        traslation[0] = -1 * ACanonical[0][3] / ACanonical[0][0];
        ACanonical[3][3] = A[3][3] - (pow(ACanonical[0][3], 2) / ACanonical[0][0]);
        ACanonical[0][3] = ACanonical[3][0] = 0;
        matrixByVector(isometry,traslation,3);/*rotating the traslation, multiplying the rotation (before multiplied by the A matrix to diagonalize the quadratic part of it) by the vector of traslation*/

        /*updating the isometry matrix with the traslation and initializating traslation to 0*/
        for(i = 0; i < 3; i++)
        {
            isometry[i][3] += traslation[i];
            traslation[i] = 0;
        }

        if(ACanonical[1][3] != 0 )
        {
            /*finding and applying the rotation that zeroes the first grade term in y*/
            double sin,cos,t, temp; //(t = sin/cos)
            if(ACanonical[2][3] == 0)
            {
                sin = 1;
                cos = 0;
            }
            else
            {
                t = ACanonical[1][3] / (-1 * ACanonical[2][3]); /*tan of the angle of the rotation needed to zero the first grade term in y*/
                cos = 1 / sqrt(pow(t, 2) + 1);
                sin = t * cos;
            }
            /*the matrix that represents this rotation is:
             * |1  0   0  |
             * |0 cos -sin|
             * |0 sin  cos|
             * */

            /*multiplying the rotation found with the jacobi method (and stored in the isometry matrix) by the new rotation that zeroes the first grade term in y */
            for(i = 0; i < 3; i++)
            {
                temp = isometry[i][1];
                isometry[i][1] = isometry[i][1] * cos + sin * isometry[i][2];
                isometry[i][2] = -sin * temp  + cos * isometry[i][2];
                
            }
            
            /*rotating the ACanonical*/
            ACanonical[2][3] = ACanonical[3][2] = -ACanonical[1][3] * sin + ACanonical[2][3] * cos;
            ACanonical[1][3] = ACanonical[3][1] = 0;

        }
        if(ACanonical[2][3] != 0)
        {
            /*traslating the quadric to zero the ACanonical[3][3]*/
            traslation[2]= -1*ACanonical[3][3] / (2 * ACanonical[2][3]);

            matrixByVector(isometry, traslation, 3);

            for(i = 0; i < 3; i++)
            {
                isometry[i][3] += traslation[i];
            }

            ACanonical[3][3] = 0;
        }
        isometry[3][3] = 1;
                
        
        
    }
}



