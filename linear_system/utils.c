#include "math.h"
#include "stdio.h"

float* getel(float* A, int n, int row, int col) {
    return A + row * n + col;
}

void print_solution(float* x, int n, int m) {
    printf("Solution: ");
    for (int i = 0; i < n && i < m; ++i) {
        printf("%f ", x[i]);
    }
    printf("\n");
}


void print_error(float* A, int n, float* b, float* x) {
    float error = 0;
    for (int row = 0; row < n; ++row) {
        float my_b = 0;
        for (int col = 0; col < n; ++col) {
            my_b += *getel(A, n, row, col) * x[col];
        }
        error += (my_b - b[row]) * (my_b - b[row]);
    }
    printf("|Ax-b|=%f\n", sqrt(error));
}


void print_system(float* A, float* b, int n) {
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            printf("%f ", *getel(A, n, row, col));
        }
        printf("| %f\n", b[row]);
    }
}


void print_solution_error(float* x, int n) {
    float solution_error = 0;
    for (int i = 0; i < n; ++i) {
        solution_error += ((1 - (i % 2)) - x[i]) * ((1 - (i % 2)) - x[i]);
    }
    printf("Solution error %f\n", sqrt(solution_error));
}


int is_zero(float x) {
    return fabs(x) < 1e-8;
}