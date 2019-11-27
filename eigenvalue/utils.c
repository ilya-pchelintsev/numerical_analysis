#include "math.h"
#include "stdio.h"

#include "utils.h"

double* getel(double* A, int n, int row, int col) {
    if (col >= row)
        return A + row * n + col;
    return A + col * n + row;
}

void print_solution(double* x, int n, int m) {
    printf("Solution: ");
    for (int i = 0; i < n && i < m; ++i) {
        printf("%f ", x[i]);
    }
    printf("\n");
}


void print_error(double* A, int n, double* b, double* x) {
    double error = 0;
    for (int row = 0; row < n; ++row) {
        double my_b = 0;
        for (int col = 0; col < n; ++col) {
            my_b += *getel(A, n, row, col) * x[col];
        }
        error += (my_b - b[row]) * (my_b - b[row]);
    }
    printf("|Ax-b|=%e\n", sqrt(error));
}


void print_matr(double* A, int n) {
    printf("np.array([");
    for (int row = 0; row < n; ++row) {
        printf("[");
        for (int col = 0; col < n - 1; ++col) {
            printf("%.10f\t, ", *getel(A, n, row, col));
        }
        printf("%.10f],\\\n", *getel(A, n, row, n - 1));
    }
    printf("])\n");
}


void print_solution_error(double* x, int n) {
    double solution_error = 0;
    for (int i = 0; i < n; ++i) {
        solution_error += ((1 - (i % 2)) - x[i]) * ((1 - (i % 2)) - x[i]);
    }
    printf("Solution error %e\n", sqrt(solution_error));
}


int is_zero(double x, double A_norm) {
    return fabs(x) < 1e-10 * A_norm;
}


double matr_norm(double* A, int n) {
    double norm = 0;
    for (int i = 0; i < n * n; ++i) {
        norm += A[i] * A[i];
    }
    return sqrt(norm);
}
