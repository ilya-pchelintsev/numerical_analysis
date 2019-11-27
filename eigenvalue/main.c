#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#include "linsolve.h"
#include "matrread.h"
#include "utils.h"

#define MAX_PRINT_SIZE 10

int main(int argc, char** argv) {
    int n, num_eigenvalues;
    double interval_beg, interval_end, precision;
    double *A, *b, *x, *diag;
    clock_t start, end;

    if (argc <= 4 || argc > 6) {
        printf("Wrong command line arguments\n");
        return 0;
    }

    if (sscanf(argv[1], "%d", &n) == 0) {
        printf("First argument should be matrix size\n");
        return 0;
    }

    if (sscanf(argv[argc - 3], "%lf", &interval_beg) == 0 ||
        sscanf(argv[argc - 2], "%lf", &interval_end) == 0 ||
        sscanf(argv[argc - 1], "%lf", &precision) == 0 ||
        precision <= 0 || interval_end < interval_beg) {
        printf("%lf\n", precision);
        printf("Wrong arguments\n");
        return 0;
    }

    A = (double*)malloc(n * n * sizeof(double));
    if (!A) {
        printf("Can't allocate memory\n");
        return 0;
    }


    if (argc == 4) {
        read_matr_formula(A, n);
    } else if (read_matr_file(A, n, argv[2]) != 0){
        printf("Can't read matrix from file");
        free(A);
        return 0;
    }

    x = (double*)malloc(n * sizeof(double));
    start = clock();
    diag = (double*)malloc(n * sizeof(double));
    num_eigenvalues = find_eigenvalues(A, n, x, interval_beg, interval_end, precision, diag);
    if (num_eigenvalues < 0) {
        printf("Can't find eigenvalues\n");
        free(A);
        free(x);
        free(diag);
        return 0;
    }
    end = clock();

    printf("Found %d eigenvalues in %f miliseconds\n",
           num_eigenvalues, ((double) (end - start)) / CLOCKS_PER_SEC * 1000);
    print_solution(x, num_eigenvalues, MAX_PRINT_SIZE);
    
    free(A);
    free(b);
    free(diag);
    return 0;
}