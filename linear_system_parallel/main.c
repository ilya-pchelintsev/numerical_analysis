#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "pthread.h"

#include "linsolve.h"
#include "matrread.h"
#include "utils.h"

#define MAX_PRINT_SIZE 10

int main(int argc, char** argv) {
    int n, num_threads;
    double *A, *b, *x, *sins, *coss;
    pthread_t *threads;
    double start_time, end_time;
    struct linsolve_args *args;

    if (argc <= 2 || argc > 4) {
        printf("Wrong command line arguments\n");
        return 0;
    }

    if (sscanf(argv[1], "%d", &n) == 0 || n <= 0) {
        printf("First argument should be matrix size\n");
        return 0;
    }

    if (sscanf(argv[argc - 1], "%d", &num_threads) == 0 ||
        num_threads <= 0) {
        printf("Last argument should be number of threads\n");
        return 0;
    }

    A = (double*)malloc(n * n * sizeof(double));
    if (!A) {
        printf("Can't allocate memory\n");
        return 0;
    }
    b = (double*)malloc(n * sizeof(double));
    if (!b) {
        printf("Can't allocate memory\n");
    }

    if (argc == 3) {
        read_matr_formula(A, b, n);
    } else if (read_matr_file(A, b, n, argv[2]) != 0){
        printf("Can't read matrix from file");
        free(A);
        free(b);
        return 0;
    }

    threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    x = (double*)malloc(n * sizeof(double));
    args = (struct linsolve_args*)malloc(num_threads * sizeof(struct linsolve_args));
    sins = (double*)malloc(n * sizeof(double));
    coss = (double*)malloc(n * sizeof(double));

    start_time = get_time();
    for (int i = 0; i < num_threads; i++) {
        args[i].A = A;
        args[i].b = b;
        args[i].n = n;
        args[i].thread_num = i;
        args[i].total_thread_num = num_threads;
        args[i].x = x;
        args[i].sins = sins;
        args[i].coss = coss;

        if (pthread_create(threads + i, 0, solve_linear_system_parallel, (void*)(args + i))) {
            printf("Can't create thread\n");
            free(A);
            free(x);
            free(b);
            free(args);
            free(threads);
            free(coss);
            free(sins);
            return 0;
        }
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], 0);
    }

    if (args[0].status != 0) {
       printf("Can't solve linear system\n");
       free(A);
       free(b);
       free(x);
       free(coss);
       free(sins);
       return 0;
    }
    end_time = get_time();

    printf("Solved linear system in %f seconds\n", end_time - start_time);
    print_solution(x, n, MAX_PRINT_SIZE);
    
    if (argc == 3) {
        print_solution_error(x, n);
    }

    if (argc == 3) {
        read_matr_formula(A, b, n);
    } else if (read_matr_file(A, b, n, argv[2]) != 0){
        printf("Can't read matrix from file to compute |Ax-y|");
        free(A);
        free(b);
        free(coss);
        free(sins);
        return 0;
    }
    print_error(A, n, b, x);

    free(A);
    free(x);
    free(args);
    free(threads);
    free(b);
    free(coss);
    free(sins);
    return 0;
}
