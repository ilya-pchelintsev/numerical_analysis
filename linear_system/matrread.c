#include "utils.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "matrread.h"

void set_b_odd_sum(double* A, double* b, int n);
double f(int i, int j);

void set_b_odd_sum(double* A, double* b, int n) {
    for (int row = 0; row < n; ++row) {
        b[row] = 0;
        for (int col = 0; col < n; col += 2) {
            b[row] += *getel(A, n, row, col);
        }
    }
}


int read_matr_file(double* A, double* b, int n, char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        return 1;
    }

    for (int i = 0; i < n * n; ++i) {
        if (fscanf(file, "%lf", A + i) != 1) {
            fclose(file);
            return 2;
        }
    }
    fclose(file);

    set_b_odd_sum(A, b, n);
    return 0;
}

double f(int i, int j) {
    return abs(i - j);
}

int read_matr_formula(double* A, double* b, int n) {
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            *getel(A, n, row, col) = f(row, col);
        }
    }

    set_b_odd_sum(A, b, n);
    return 0;
}
