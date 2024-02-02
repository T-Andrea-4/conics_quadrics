#include <stdio.h>
#include <stdlib.h>
#include "gauss.h"
#include "conic.h"

int main() {
    double **a = (double**) malloc(sizeof(double *) * 3);
    *a = (double *) malloc(sizeof(double ) * 3 * 3);
    for(int i = 0; i < 3; i++ )
    {
        *(a + i) = *(a) + 3 * i;
    }

    a[0][0] = 25;
    a[0][1] = 0;
    a[0][2] = 0;
    a[1][0] = 0;
    a[1][1] = -7;
    a[1][2] = 24;
    a[2][0] = 0;
    a[2][1] = 24;
    a[2][2] = 7;

    double **aCanonical = (double**) malloc(sizeof(double *) * 3);
    *aCanonical = (double *) malloc(sizeof(double ) * 3 * 3);
    for(int i = 0; i < 3; i++ )
    {
        *(aCanonical + i) = *(aCanonical) + 3 * i;
    }

    double **isometria = (double**) malloc(sizeof(double *) * 3);
    *isometria = (double *) malloc(sizeof(double ) * 3 * 3);
    for(int i = 0; i < 3; i++ )
    {
        *(isometria + i) = *(isometria) + 3 * i;
    }

    double *centre = (double*) malloc(sizeof (double) * 2);

    canonicalForm(a,aCanonical, isometria, centre);

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf(" %f", aCanonical[i][j]);
        }
        printf("\n");
    }

    printf("\n%f, %f", centre[0], centre[1]);

    return 0;

}