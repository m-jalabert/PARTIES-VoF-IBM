/* L2 norm */
#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i;
    int LEN = 100; // node points in y direction
    float *v1 = (float *)malloc(LEN * sizeof(float));
    float *v2 = (float *)malloc(LEN * sizeof(float));
    float *result = (float *)malloc(LEN * sizeof(float));
    char c[10];
    FILE *fptr;

    fptr = fopen("analytical_velo.dat", "r");
    if (fptr)
    {
        // Reading analytical velocities
        for (i = 0; i < LEN; i++)
        {
            if (i == 0)
                fscanf(fptr, "%s", c);

            fscanf(fptr, "%f", &v1[i]);
        }
        // Print analytical velocities
        printf("Elements of analytical velocity:\n");
        for (i = 0; i < LEN; i++)
            printf("%f\n", v1[i]);
        printf("\n");
    }
    else
    {
        printf("Error! opening file");
        // Exception: Program exits if file pointer returns NULL.
    }

    fptr = fopen("numerical_velo.dat", "r");
    if (fptr)
    {
        // Reading numerical velocities
        for (i = 0; i < LEN; i++)
        {
            if (i == 0)
                fscanf(fptr, "%s", c);

            fscanf(fptr, "%f", &v2[i]);
        }

        // Print numerical velocities
        printf("Elements of numerical velocity:\n");
        for (i = 0; i < LEN; i++)
            printf("%f\n", v2[i]);
        printf("\n");
    }
    else
    {
        printf("Error! opening file");
        // Exception: Program exits if file pointer returns NULL.
    }

    // l2 norm
    for (i = 0; i < LEN; i++)
        result[i] = 0.5 * ((v1[i] - v2[i]) * (v1[i] - v2[i]));

    // Print result vector
    printf("Resultant vector:\n");
    for (i = 0; i < LEN; i++)
        printf("%f\n", result[i]);
    printf("\n");

    fclose(fptr);
    free(v1);
    free(v2);
    free(result);
    return 0;
}