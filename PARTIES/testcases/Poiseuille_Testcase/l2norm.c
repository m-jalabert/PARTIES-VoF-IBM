/* C program to calculate the L2 norm 
* between two diffect vectors present in .dat files */
#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i=0;
    int testresult = 1; // true
    int LEN = 101; // node points in y direction
    float *y = (float *)malloc(LEN * sizeof(float));
    float *v1 = (float *)malloc(LEN * sizeof(float));
    float *v2 = (float *)malloc(LEN * sizeof(float));
    float *result = (float *)malloc(LEN * sizeof(float));
    char c[10];
    char d[10];
    FILE *fptr;

    fptr = fopen("numerical_velo_verified.dat", "r");
    if (fptr)
    {
        while(!feof(fptr))
        {
            if (i == 0)
            {
                fscanf(fptr, "%s\t%s", c, d);
            }
            else
            {
                fscanf(fptr, "%f\t%f", &y[i-1], &v1[i-1]);
            }
            i++;
        }

        // Print analytical velocities
        // printf("Elements of verified velocity:\n");
        // for (i = 0; i < LEN; i++)
        //     printf("%10.8f\n", v1[i]);
        // printf("\n");
    }
    else
    {
        printf("Error! opening file");
        return 1;
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
        // printf("Elements of numerical velocity:\n");
        // for (i = 0; i < LEN; i++)
        //     printf("%f\n", v2[i]);
        // printf("\n");
    }
    else
    {
        printf("Error! opening file");
        return 1;
        // Exception: Program exits if file pointer returns NULL.
    }

    // l2 norm
    for (i = 0; i < LEN; i++)
        result[i] = 0.5 * ((v1[i] - v2[i]) * (v1[i] - v2[i]));

    // writing l2nor to file
    fptr = fopen("l2norm.dat", "w");
    if (fptr)
    {
        // Save l2norm
        for (i = 0; i < LEN; i++)
        {
            if (i == 0)
                fprintf(fptr, "%s\n", "      U(y)");

            fprintf(fptr, "%f\n", result[i]);
        }
    }
    else
    {
        printf("Error! opening file");
        return 1;
        // Exception: Program exits if file pointer returns NULL.
    }

    // Print result vector
    // printf("Resultant L2 norm of vectors:\n");
    // for (i = 0; i < LEN; i++)
    //     printf("%f\n", result[i]);
    // printf("\n");

    // Check if pass or fail
    for (i = 0; i < LEN; i++)
    {
        if (result[i] < 0.000000009)
        {
            continue;
        }
        else
        {
            testresult = 0;
            break;
        }
    }

    if (testresult == 1)
    {
        printf("----------**Test PASS!!**----------\n");
        fprintf(fptr, "%s\n", "PASS");
    }
    else
    {
        printf("----------**Test FAIL!!**----------\n");
        fprintf(fptr, "%s\n", "FAIL");
    }

    // Adding coordinates to numerical_velo.dat
    fptr = fopen("numerical_velo.dat", "w");
    if (fptr)
    {
        // Save cood and numerical
        for (i = 0; i <= LEN; i++)
        {
            if (i == 0)
                fprintf(fptr, "%s\t%s\n", "y","U(y)");
            else
            {
                fprintf(fptr, "%5.4f\t%10.8f\n", y[i-1],v2[i-1]);
            }
        }
    }
    else
    {
        printf("Error! opening file");
        return 1;
        // Exception: Program exits if file pointer returns NULL.
    }

    fclose(fptr);
    free(y);
    free(v1);
    free(v2);
    free(result);
    return 0;
}
