#include <stdio.h>
#include "util.h"

void print_array(const int A[], const int size)
{
    int i = 0;

    while (i < size)
        printf("%d ", A[i++]);
    printf("\n");
}

void print_float_array(const float A[], const int size)
{
    int i = 0;

    while (i < size)
        printf("%.2f ", A[i++]);
    printf("\n");
}

void print_array_to_file(FILE *fp, const int A[], const int size)
{
    int i = 0;

    while (i < size)
        fprintf(fp, "%d ", A[i++]);
    fprintf(fp, "\n");
}

void print_float_array_to_file(FILE *fp, const float A[], const int size)
{
    int i = 0;

    while (i < size)
        fprintf(fp, "%.2f ", A[i++]);
    fprintf(fp, "\n");
}
