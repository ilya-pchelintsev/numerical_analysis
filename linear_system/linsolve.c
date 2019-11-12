#include "utils.h"
#include "stdio.h"
#include "math.h"

void turn_lines(float* A, float* b, int n, int row1, int row2, int start_col) {
    float c = *getel(A, n, row1, start_col);
    float s = *getel(A, n, row2, start_col);

    if (is_zero(s)) {
        return;
    }

    float norm = sqrt(c * c + s * s);
    c /= norm;
    s /= norm;

    float* cur_row1 = getel(A, n, row1, start_col);
    float* cur_row2 = getel(A, n, row2, start_col);
    float old_row1, old_row2;
    for (int col = start_col; col < n; ++col) {        
        old_row1 = *cur_row1;
        old_row2 = *cur_row2;
        *cur_row1 =  c * old_row1 + s * old_row2;
        *cur_row2 = -s * old_row1 + c * old_row2;
        ++cur_row1;
        ++cur_row2;
    }

    float old_b1 = b[row1];
    float old_b2 = b[row2];
    b[row1] =  c * old_b1 + s * old_b2;
    b[row2] = -s * old_b1 + c * old_b2;
}


void triangular_form(float* A, float* b, int n) {
    for (int col = 0; col < n; ++col) {
        for (int row = col + 1; row < n; ++row) {
            turn_lines(A, b, n, col, row, col);
        }
    }
}


void add_line(float* A, float* b, int n, int dest_row, int src_row, float mult, int start_col) {
    float* src_row_ptr = getel(A, n, src_row, start_col);
    float* dest_row_ptr = getel(A, n, dest_row, start_col);
    for (int col = start_col; col < n; ++col) {
        *dest_row_ptr += mult * (*src_row_ptr);
        ++src_row_ptr;
        ++dest_row_ptr;
    }

    b[dest_row] += mult * b[src_row];
}


int diagonal_form(float* A, float* b, int n) {
    for (int col = n - 1; col >= 0; --col) {
        float diag_el = *getel(A, n, col, col);
        if (is_zero(diag_el)) {
            return 1;
        }
        for (int row = col - 1; row >= 0; --row) {
            float mult = -(*getel(A, n, row, col)) / diag_el;
            add_line(A, b, n, row, col, mult, col);
        }
    }
    return 0;
}


int solve_linear_system(float* A, int n, float* b, float* x) {
    triangular_form(A, b, n);

    if (diagonal_form(A, b, n) != 0) {
        return 1;
    }

    for (int i = 0; i < n; ++i) {
        x[i] = b[i] / (*getel(A, n, i, i));
    }
    return 0;
}