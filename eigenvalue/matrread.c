#include "utils.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "matrread.h"

double f(int i, int j);


int read_matr_file(double* A, int n, char* filename) {
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
    return 0;
}

double f(int i, int j) {
    return abs(i - j);
}

int read_matr_formula(double* A, int n) {
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            *getel(A, n, row, col) = f(row, col);
        }
    }
    return 0;
}
