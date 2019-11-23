#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#include "linsolve.h"
#include "matrread.h"
#include "utils.h"

#define MAX_PRINT_SIZE 10

int main(int argc, char** argv) {
    int n;
    double *A, *b, *x;
    clock_t start, end;

    if (argc <= 1 || argc > 3) {
        printf("Wrong command line arguments\n");
        return 0;
    }

    if (sscanf(argv[1], "%d", &n) == 0) {
        printf("First argument should be matrix size\n");
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

    if (argc == 2) {
        read_matr_formula(A, b, n);
    } else if (read_matr_file(A, b, n, argv[2]) != 0){
        printf("Can't read matrix from file");
        free(A);
        free(b);
        return 0;
    }

    x = (double*)malloc(n * sizeof(double));
    start = clock();
    if (solve_linear_system(A, n, b, x) != 0) {
        printf("Can't solve linear system\n");
        free(A);
        free(b);
        free(x);
        return 0;
    }
    end = clock();

    printf("Solved linear system in %f miliseconds\n",
           ((double) (end - start)) / CLOCKS_PER_SEC * 1000);
    print_solution(x, n, MAX_PRINT_SIZE);
    
    if (argc == 2) {
        print_solution_error(x, n);
    }

    if (argc == 2) {
        read_matr_formula(A, b, n);
    } else if (read_matr_file(A, b, n, argv[2]) != 0){
        printf("Can't read matrix from file to compute |Ax-y|");
        free(A);
        free(b);
        return 0;
    }
    print_error(A, n, b, x);

    free(A);
    free(x);
    free(b);
    return 0;
}