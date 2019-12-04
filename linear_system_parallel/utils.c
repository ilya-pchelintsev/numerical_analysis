#include "math.h"
#include "stdio.h"
#include "pthread.h"
#include "sys/time.h"

#include "utils.h"

double* getel(double* A, int n, int row, int col) {
    return A + row * n + col;
}

void print_solution(double* x, int n, int m) {
    printf("Solution: ");
    for (int i = 0; i < n && i < m; ++i) {
        printf("%.4f ", x[i]);
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


void print_system(double* A, double* b, int n) {
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            printf("%f ", *getel(A, n, row, col));
        }
        printf("| %f\n", b[row]);
    }
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

void synchronize(int total_threads)
{
    //printf("sync start\n");
    static pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in=PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out=PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    pthread_mutex_lock (&mutex);
    threads_in++;
    if (threads_in >= total_threads) {
        threads_out = 0;
        pthread_cond_broadcast (&condvar_in);
    }
    else {
        while (threads_in < total_threads) {
            pthread_cond_wait (&condvar_in, &mutex);
        }
    }
    threads_out++;
    if (threads_out >= total_threads) {
        threads_in = 0;
        pthread_cond_broadcast (&condvar_out);
    }
    else {
        while (threads_out < total_threads) {
            pthread_cond_wait (&condvar_out, &mutex);
        }
    }
    pthread_mutex_unlock (&mutex);
    //printf("sync end\n");
}

double get_time(void) {
    struct timeval time;
    gettimeofday(&time, 0);
    return (double)(time.tv_sec + time.tv_usec / 1000000.0);
}
