#include "utils.h"
#include "stdio.h"
#include "math.h"

#include "linsolve.h"

void turn_lines(double* A, double* b, int n, int row1, int row2, int start_col) {
    double *cur_row1, *cur_row2, *end_row1;
    double old_row1, old_row2, norm, c, s, old_b1, old_b2;

    c = *getel(A, n, row1, start_col);
    s = *getel(A, n, row2, start_col);

    if (is_zero(s)) {
        return;
    }

    norm = sqrt(c * c + s * s);
    c /= norm;
    s /= norm;

    cur_row1 = getel(A, n, row1, start_col);
    end_row1 = getel(A, n, row1 + 1, 0);
    cur_row2 = getel(A, n, row2, start_col);
    while (cur_row1 < end_row1) {
        old_row1 = *cur_row1;
        old_row2 = *cur_row2;
        *cur_row1 =  c * old_row1 + s * old_row2;
        *cur_row2 = -s * old_row1 + c * old_row2;
        ++cur_row1;
        ++cur_row2;
    }

    old_b1 = b[row1];
    old_b2 = b[row2];
    b[row1] =  c * old_b1 + s * old_b2;
    b[row2] = -s * old_b1 + c * old_b2;
}


void triangular_form(double* A, double* b, int n) {
    for (int col = 0; col < n; ++col) {
        for (int row = col + 1; row < n; ++row) {
            turn_lines(A, b, n, col, row, col);
        }
    }
}


void add_line(double* A, double* b, int n, int dest_row, int src_row, double mult, int start_col) {
    b[dest_row] += mult * b[src_row];
}


int diagonal_form(double* A, double* b, int n) {
    for (int col = n - 1; col >= 0; --col) {
        double diag_el = *getel(A, n, col, col);
        if (is_zero(diag_el)) {
            return 1;
        }
        for (int row = col - 1; row >= 0; --row) {
            double mult = -(*getel(A, n, row, col)) / diag_el;
            add_line(A, b, n, row, col, mult, col);
        }
    }
    return 0;
}


int solve_linear_system(double* A, int n, double* b, double* x) {
    triangular_form(A, b, n);

    if (diagonal_form(A, b, n) != 0) {
        return 1;
    }

    for (int i = 0; i < n; ++i) {
        x[i] = b[i] / (*getel(A, n, i, i));
    }
    return 0;
}
