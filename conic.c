#include <math.h>
#include "gauss.h"
#include "diagonalization.h"


/* in this file B is the matrix that represents the quadratic part (quadratic form) of the conic section while
 * A is the complete 3 X 3 matrix associated to the conic section. To find the canonical form of the conic,
 * the algorithm uses the othogonal invariants method to recognize the conic section*/


char recognition(double B[2][2])
{
    /*the function returns 'p','h','e' if the conic shape is respectively a parabola,
     * a hyperbola, an ellipse.*/
    double detB;
    detB = determinant(2, 2,B);
    if(detB > 0)
    {
        return 'e';
    }
    else if(detB == 0)
    {
        return 'p';
    }
    else
    {
        return 'h';
    }


}

void conicCanonicalForm(double A[3][3], double ACanonical[3][3], double isometry[3][3], double centreOrVertex[2], char *type)
{
    /*If the conic is degenerate  'd is returned and others
     * arguments are not modified.*/
    if(determinant(3,3,A) == 0)
    {
        *type = 'd';
        return;
    }
    
    /*isometry is the matrix that represents the translation and the rotation that transform the
     * conic shape in the canonical form*/

    int i, j;
    double B[2][2];
    double rotation[2][2];
    
    /*initializing B as the matrix of the quadratic part of the conic*/
    for( i = 0; i < 2; i++)
    {
        for(j = 0; j < 2; j++)
        {
            B[i][j] = A[i][j];
        }
    }
    
    /*recognizing conic type*/
    *type = recognition(B);

    
    
    if(*type != 'p') //the conic has a centre
    {
        double centreSystem[2][3]; //this 2 X 3 matrix represents the system that must be solved to find the centre
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 3; j++)
            {
                if(j == 2)
                {
                    centreSystem[i][j] = -1 * A[i][j];
                }
                else{
                    centreSystem[i][j] = A[i][j];
                }
            }
        }
        /*gauss and backword will solve the system*/
        gauss(2, 3, centreSystem);
        backward(2,3, centreSystem);
        centreOrVertex[0] = centreSystem[0][2];
        centreOrVertex[1] = centreSystem[1][2];

        jacobi(2,B,ACanonical,rotation); //finds the rotation that diagonalizes the quadratic part

        ACanonical[2][2] = determinant(3,3,A) / (ACanonical[0][0] * ACanonical[1][1]); //initialization of the constant term
        ACanonical[0][2] = ACanonical[2][0] = ACanonical[1][2] = ACanonical[2][1] = 0; // setting all the elements out of the diagonal to 0
    }
    else
    {
        double v[2];/*this will be set as the vertex of the parabola (v[0] = x, v[1] = y) after the rotation (will be used to find the real vertex not the one after the rotation that will be (0,0)*/
        double a,b,c;

        jacobi(2, B, ACanonical, rotation);//finds the rotation that diagonalizes the quadratic part

        /*changing the order of the elements in the diagonal of B if ACanonical[0][0] = 0 (and consequently changing the order of the eigenvectors in the rotation matrix) this is done so that the parabla has the
         y =aX^2 + bX + c*/
        if(ACanonical[0][0] == 0)
        {
            ACanonical[0][0] = ACanonical[1][1];
            ACanonical[1][1] = 0;
            columnSwap(2,2,rotation,0,1);
        }

        /*applying the rotation found with the jacobi method to the linear part of the canonical matrix (that then will be used to find the vertex of the parabola*/
        ACanonical[0][2] = ACanonical[2][0] = rotation[0][0] * A[0][2] + rotation[0][1] * A[1][2];
        ACanonical[1][2] = ACanonical[2][1] = rotation[1][0] * A[0][2] + rotation[1][1] * A[1][2];
        ACanonical[2][2] = A[2][2];

        /*given y = aX^2 + bX + c, the x coordinate of the vertex is -b/2a*/
        a = ACanonical[0][0] / (-2 * ACanonical[1][2]);
        b = -1 * ACanonical[0][2] / ACanonical[1][2];
        c = -0.5 * (ACanonical[2][2] / ACanonical[1][2]);

        v[0] = -1 * b / (2*a);
        v[1] = a * pow(v[0],2) + b * v[0] + c;

        /* finding the vertex of the parabola in the original frame of reference (using the inverse eigenvectors matrix) */
        centreOrVertex[0] = rotation[0][0] * v[0] + rotation[1][0] * v[1];
        centreOrVertex[1] = rotation[0][1] * v[0] + rotation[1][1] * v[1];
        ACanonical[0][2] = ACanonical[2][0] = 0;
        ACanonical[2][2] = 0;
        ACanonical[1][2] = ACanonical[2][1] = -1 *sqrt(-1 * (determinant(3, 3, A)) / ACanonical[0][0]);
    }

    isometry[2][0] = isometry[2][1] = 0;
    /*translation part of the isometry*/
    isometry[0][2] = centreOrVertex[0];
    isometry[1][2] = centreOrVertex[1];
    isometry[2][2] = 1;
    
    /*rotation part of the isometry*/
    for(i = 0; i < 2; i++ )
    {
        for(j = 0; j < 2; j++)
        {
            isometry[i][j] = rotation[i][j];
        }
    }
}
