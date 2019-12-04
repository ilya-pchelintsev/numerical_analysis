#include "pthread.h"

struct linsolve_args {
    double *A;
    int n;
    double *b;
    double* x;
    int thread_num;
    int total_thread_num;
    int status;
    //struct pthread_barrier_t* barrier;
};

void solve_linear_system_parallel(void *args);

//int solve_linear_system(double* A, int n, double* b, double* x);