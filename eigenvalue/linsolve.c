#include "utils.h"
#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#include "linsolve.h"

#define EPS 1e-15


static void turn_lines_rows(double* A, int n, int row1, int row2, double A_norm) {
    double *cur_row1, *cur_row2, *end_row1;
    double old_row1, old_row2, norm, c, s, old_b1, old_b2, x, y, z;

    c = *getel(A, n, row1, row1 + 1);
    s = *getel(A, n, row1, row2);

    if (fabs(s * A_norm) < EPS) {
        return;
    }

    norm = sqrt(c * c + s * s);
    c /= norm;
    s /= norm;

    for (int cur_row = row1; cur_row < n; cur_row++) {
        if (cur_row == row1 + 1 || cur_row == row2) {
            continue;
        }
        cur_row1 = getel(A, n, cur_row, row1 + 1);
        cur_row2 = getel(A, n, cur_row, row2);
        old_row1 = *cur_row1;
        old_row2 = *cur_row2;
        *cur_row1 =  c * old_row1 + s * old_row2;
        *cur_row2 = -s * old_row1 + c * old_row2;
    }

    x = *getel(A, n, row1 + 1, row1 + 1);
    y = *getel(A, n, row1 + 1, row2);
    z = *getel(A, n, row2, row2);

    *getel(A, n, row1 + 1, row1 + 1) = c * c * x + 2 * s * c * y + s * s * z;
    *getel(A, n, row1 + 1, row2) = - s * c * x + s * c * z + c * c * y - s * s * y;
    *getel(A, n, row2, row2) = s * s * x + -2 * s * c * y + c * c * z;

    //printf("%lf %lf\n", c, s);
}

int lu_diag(double* A, int n, double* diag, double A_norm) {
    double u;

    diag[0] = *getel(A, n, 0, 0);

    if (fabs(diag[0]) < A_norm * EPS) {
        return -1;
    }

    u = *getel(A, n, 0, 1) / diag[0];
    for (int i = 1; i < n; i++) {
        diag[i] = *getel(A, n, i, i) - *getel(A, n, i, i - 1) * u;
        if (fabs(diag[i]) < A_norm * EPS) {
            return -1;
        }
        if (i < n - 1)
            u = *getel(A, n, i + 1, i) / diag[i];
    }
    return 0;
}

void tridiagonal_form(double* A, int n, double A_norm) {
    //print_matr(A, n);
    //printf("----------------------------------\n");
    for (int i = 0; i < n - 2; i++) {
        for (int j = i + 2; j < n; j++) {
            turn_lines_rows(A, n, i, j, A_norm);
            //print_matr(A, n);
            //printf("----------------------------------\n");
        }
    }

    //print_matr(A, n);
    //printf("----------------------------------\n");
}

int num_sign_changes(double* x, int n, double A_norm) {
    if (fabs(x[0]) < A_norm * EPS) {
        printf("ALARM 2\n");
        return -1;
    }
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(x[i]) < A_norm * EPS) {
            printf("ALARM %d\n", i);
            return -1;
        } if (x[i] < 0)
            count++;
    }
    return count;
}

int num_eigenvalues_less(double* A, int n, double lambda, double A_norm, double* diag) {
    for (int i = 0; i < n; i++) {
        *getel(A, n, i, i) -= lambda;
    }

    if (lu_diag(A, n, diag, A_norm)){
        return -1;
    }

    //for (int i = 0; i < n; i++) {
    //    printf("%lf ", diag[i]);
    //}
    //printf("op\n");

    int count = num_sign_changes(diag, n, A_norm);
    //printf("num neg %d\n", count);

    for (int i = 0; i < n; i++) {
        *getel(A, n, i, i) += lambda;
    }

    return count;
}

int find_eigenvalues(double* A, int n, double* x, double low, double high, double precision, double* diag) {
    int low_n, high_n;
    int num_eigenvalues = 0;
    
    double A_norm = matr_norm(A, n);
    tridiagonal_form(A, n, A_norm);

    //double* diag = (double*)malloc(n * sizeof(double));
    
    low_n = num_eigenvalues_less(A, n, low, A_norm, diag);
    high_n = num_eigenvalues_less(A, n, high, A_norm, diag);

    //printf("%d %d\n", low_n, high_n);

    for (int k = low_n + 1; k <= high_n; k++) {
        double a, b;
        int n_a, n_b;
        a = low;
        b = high;
        n_a = num_eigenvalues_less(A, n, a, A_norm, diag);
        n_b = num_eigenvalues_less(A, n, b, A_norm, diag);
        while (fabs(b - a) > precision) {
            double c = (a + b) / 2;
            int n_c = num_eigenvalues_less(A, n, c, A_norm, diag);
            if (n_c < 0) {
                a -= precision;
                //printf("ALARM3\n");
                continue;
            }
            if (n_c >= k) {
                b = c;
            } else {
                a = c;
            }
            n_a = num_eigenvalues_less(A, n, a, A_norm, diag);
            n_b = num_eigenvalues_less(A, n, b, A_norm, diag);

            //printf("%d: %e %e %d %d\n", k, a, b, n_a, n_b);
        }

        //printf("found %lf mult %d\n", (a + b) / 2, n_b - n_a);
        for (int i = 0; i < n_b - n_a; i++) {
            x[num_eigenvalues] = (a + b) / 2;
            num_eigenvalues++;
        }
        k += n_b - n_a - 1;
    }

    free(diag);

    return high_n - low_n;
}
