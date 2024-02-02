#include <stdlib.h>
#include <math.h>
#include "gauss.h"

/* in this file B is the matrix that represents the quadratic part of the conic section while
 * A is the matrix associated to the conic section. To find the canonical form of the conic,
 * the algorithm uses the othogonal invariants method*/

void jacobi(double **B, double **ACanonical, double **isometry);


char recognition(double **A)
{
    /*the function returns 'p','h','e' if the conic shape is respectively a parabola,
     * a hyperbola, an ellipse*/
    double detB;
    detB = determinant(A, 2,2);
    if(detB > 0)
    {
        return 'e';
    }
    else if(detB == 0.0)
    {
        return 'p';
    }
    else
    {
        return 'h';
    }
    
    
}

void canonicalForm(double **A, double **ACanonical, double **isometry, double *centre)
{
    /*isometry is the matrix that represents the translation and the rotation that bring the
     * conic shape in the canonical form*/

    char conicType;
    int i, j;
    conicType = recognition(A);

    if(conicType != 'p') //the conic has a centre
    {
        double **centreSystem = (double**) malloc(sizeof(double *) * 2);  //system that solved gives the coordinates of the centre of the conic
        *centreSystem = (double *) malloc(sizeof(double ) * 2 * 3);
        for(i = 0; i < 2; i++ )
        {
            *(centreSystem + i) = *(centreSystem) + 3 * i;
        }

        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 3; j++)
            {
                centreSystem[i][j] = A[i][j];
            }
        }


        for(i = 0; i < 3; i++)
        {
            for(j = 0; j < 3; j++)
            {
                if(i != j)
                {
                    ACanonical[i][j] = 0;
                }
            }
        }

        jacobi(A,ACanonical,isometry);

        ACanonical[2][2] = determinant(A,3,3) / (ACanonical[0][0] * ACanonical[1][1]);

        backward(centreSystem, 2);
        centre[0] = centreSystem[0][2];
        centre[1] = centreSystem[1][2];

        isometry[2][2] = 1;
        isometry[0][2] = centre[0];
        isometry[1][2] = centre[1];
        isometry[2][0] = 0;
        isometry[2][1] = 0;
    }
    else
    {

    }
}

int sign(double d)
{
    if(d > 0)
    {
        return 1;
    }
    if(d < 0)
    {
        return -1;
    }
    return 0;
}


void jacobi(double **B, double **ACanonical, double **isometry)
{
    double theta, t, c,s;
    theta = (B[1][1] - B[0][0]) / (2 * B[0][1]);
    t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));

    ACanonical[0][1] = 0;
    ACanonical[1][0] = 0;

    ACanonical[0][0] = B[0][0] - (t * B[0][1]);
    ACanonical[1][1] = B[1][1] + (t * B[0][1]);

    c = 1 / sqrt(pow(t,2) + 1);
    s = t * c;

    isometry[0][0] = isometry[1][1] = c;
    isometry[1][0] = -1 * s;
    isometry[0][1] = s;

}



