#include "utils.h"
#include "stdio.h"
#include "math.h"
#include "pthread.h"
#include "stdlib.h"

#include "linsolve.h"
#include "utils.h"
/*
static void turn_lines(double* A, double* b, int n, int row1, int row2, int start_col, double A_norm, int thread_num, int total_thread_num) {
    double *cur_row1, *cur_row2, *end_row1;
    double old_row1, old_row2, norm, c, s, old_b1, old_b2;

    c = *getel(A, n, row1, start_col);
    s = *getel(A, n, row2, start_col);

    if (is_zero(s, A_norm)) {
        return;
    }

    norm = sqrt(c * c + s * s);
    c /= norm;
    s /= norm;

    cur_row1 = getel(A, n, row1, start_col) + thread_num + 1;
    end_row1 = getel(A, n, row1 + 1, 0);
    cur_row2 = getel(A, n, row2, start_col) + thread_num + 1;
    while (cur_row1 < end_row1) {
        old_row1 = *cur_row1;
        old_row2 = *cur_row2;
        *cur_row1 =  c * old_row1 + s * old_row2;
        *cur_row2 = -s * old_row1 + c * old_row2;
        cur_row1 += total_thread_num;
        cur_row2 += total_thread_num;
    }

    synchronize(total_thread_num);

    if (thread_num == 0) {
        cur_row1 = getel(A, n, row1, start_col);
        cur_row2 = getel(A, n, row2, start_col);
        old_row1 = *cur_row1;
        old_row2 = *cur_row2;
        *cur_row1 =  c * old_row1 + s * old_row2;
        *cur_row2 = -s * old_row1 + c * old_row2;

        old_b1 = b[row1];
        old_b2 = b[row2];
        b[row1] =  c * old_b1 + s * old_b2;
        b[row2] = -s * old_b1 + c * old_b2;
    }
}


static void rotation_method(double* A, double* b, int n, double A_norm,
                            int thread_num, int total_thread_num, int col) {
    for (int row = col + 1; row < n; row++) {
        turn_lines(A, b, n, col, row, col, A_norm, thread_num, total_thread_num);
        synchronize(total_thread_num);
    }
}


*/

static void swap_lines(double* A, double* b, int n, int row1, int row2, int start) {
    double x;
    for (int i = start; i < n; i++) {
        x = A[n * row1 + i];
        A[n * row1 + i] = A[n * row2 + i];
        A[n * row2 + i] = x;
    }

    x = b[row1];
    b[row1] = b[row2];
    b[row2] = x;
}

static void add_full_line(double* A, double* b, int n, int src, int dest) {
    double mult = A[n * dest + src] / A[n * src + src];
    for (int i = src; i < n; i++) {
        A[n * dest + i] -= mult * A[n * src + i];
    }

    b[dest] -= mult * b[src];
}

static void gauss_method(double* A, double* b, int n, double A_norm,
                         int thread_num, int total_thread_num, int col) {
    if (thread_num == 0) {
        int max_row = col;
        for (int i = col + 1; i < n; i++) {
            if (A[i * n + col] > A[n * max_row + col])
                max_row = i;
        }
        if (max_row != col) {
            swap_lines(A, b, n, col, max_row, col);
        }
    }
    synchronize(total_thread_num);

    if (fabs(A[col * n + col]) < EPS * A_norm) {
        return;
    }

    for (int row = col + 1 + thread_num; row < n; row += total_thread_num) {
        add_full_line(A, b, n, col, row);
    }
}




static void triangular_form(double* A, double* b, int n, double A_norm,
                            int thread_num, int total_thread_num) {
    for (int col = 0; col < n; ++col) {
        gauss_method(A, b, n, A_norm, thread_num, total_thread_num, col);
        synchronize(total_thread_num);

        //rotation_method(A, b, n, A_norm, thread_num, total_thread_num, col);
    }
}

static int diagonal_form(double* A, double* b, int n, double A_norm) {
    for (int col = n - 1; col >= 0; --col) {
        double diag_el = *getel(A, n, col, col);
        if (is_zero(diag_el, A_norm)) {
            return 1;
        }
        for (int row = col - 1; row >= 0; --row) {
            double mult = -(*getel(A, n, row, col)) / diag_el;
            b[row] += mult * b[col];
        }
    }
    return 0;
}


static int solve_linear_system(double* A, int n, double* b, double* x, int thread_num, int total_thread_num) {
    double A_norm = matr_norm(A, n);
    synchronize(total_thread_num);

    triangular_form(A, b, n, A_norm, thread_num, total_thread_num);

    if (thread_num == 0) {
        if (diagonal_form(A, b, n, A_norm) != 0) {
            return 1;
        }

        for (int i = 0; i < n; ++i) {
            x[i] = b[i] / (*getel(A, n, i, i));
        }
        return 0;
    }
    return 0;
}

void* solve_linear_system_parallel(void *args) {
    struct linsolve_args *linsolve_args = (struct linsolve_args*)args;
    linsolve_args->status = solve_linear_system(linsolve_args->A, linsolve_args->n, linsolve_args->b,
        linsolve_args->x, linsolve_args->thread_num, linsolve_args->total_thread_num);
    return NULL;
}
