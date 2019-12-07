#include "utils.h"
#include "stdio.h"
#include "math.h"
#include "pthread.h"
#include "stdlib.h"

#include "linsolve.h"
#include "utils.h"

static void turn_column_fill_sin_cos(double* A, double* b, int n, int col, double A_norm, 
                                     int thread_num, int total_thread_num,
                                     double* sins, double* coss) {
    double* row_start_ptr, *cur_row_ptr;
    row_start_ptr = getel(A, n, col, col);
    cur_row_ptr = getel(A, n, col + 1, col);

    for (int row = col + 1; row < n; row++) {
        double old_row1, old_row2, norm;

        norm = sqrt((*cur_row_ptr) * (*cur_row_ptr) + 
                    (*row_start_ptr) * (*row_start_ptr));
        
        if (fabs(norm) < EPS * A_norm) {
            coss[row] = 1;
            sins[row] = 0;
        } else {
            coss[row] = (*row_start_ptr) / norm;
            sins[row] = (*cur_row_ptr) / norm;
        }

        old_row1 = *row_start_ptr;
        old_row2 = *cur_row_ptr;
        *row_start_ptr =  coss[row] * old_row1 + sins[row] * old_row2;
        *cur_row_ptr =   -sins[row] * old_row1 + coss[row] * old_row2;

        old_row1 = b[col];
        old_row2 = b[row];
        b[col] =  coss[row] * old_row1 + sins[row] * old_row2;
        b[row] = -sins[row] * old_row1 + coss[row] * old_row2;

        //print_system(A, b, n);
        ++cur_row_ptr;
    }
}

static void turn_column(double* A, double* b, int n, int col, int start_row, double A_norm, 
                                     int thread_num, int total_thread_num,
                                     double* sins, double* coss) {
    double* row_start_ptr, *cur_row_ptr;
    row_start_ptr = getel(A, n, start_row, col);
    cur_row_ptr = getel(A, n, start_row + 1, col);
    
    for (int row = start_row + 1; row < n; row++) {
        double old_row1, old_row2;
        old_row1 = *row_start_ptr;
        old_row2 = *cur_row_ptr;
        *row_start_ptr =  coss[row] * old_row1 + sins[row] * old_row2;
        *cur_row_ptr =   -sins[row] * old_row1 + coss[row] * old_row2;
        ++cur_row_ptr;
    }
}


static void triangular_form(double* A, double* b, int n, double A_norm,
                            int thread_num, int total_thread_num,
                            double* sins, double* coss) {
    for (int col = 0; col < n; ++col) {
        if (thread_num == 0) {
            turn_column_fill_sin_cos(A, b, n, col, A_norm, thread_num, total_thread_num, sins, coss);
        }
        synchronize(total_thread_num);

        for (int col2 = col + 1 + thread_num; col2 < n; col2 += total_thread_num) {
            turn_column(A, b, n, col2, col, A_norm, thread_num, total_thread_num, sins, coss);
        }
        synchronize(total_thread_num);
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


static int solve_linear_system(double* A, int n, double* b, double* x, 
                               int thread_num, int total_thread_num,
                               double* sins, double* coss) {
    double A_norm = matr_norm(A, n);
    synchronize(total_thread_num);

    triangular_form(A, b, n, A_norm, thread_num, total_thread_num, sins, coss);

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
        linsolve_args->x, linsolve_args->thread_num, linsolve_args->total_thread_num,
        linsolve_args->sins, linsolve_args->coss);
    return NULL;
}
