#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#include "linsolve.h"
#include "matrread.h"
#include "utils.h"

#define MAX_PRINT_SIZE 10

int main(int argc, char** argv) {
    if (argc <= 1 || argc > 3) {
        printf("Wrong command line arguments\n");
        return 0;
    }

    int n;
    if (sscanf(argv[1], "%d", &n) == 0) {
        printf("First argument should be matrix size\n");
        return 0;
    }

    float* A = (float*)malloc(n * n * sizeof(float));
    float* b = (float*)malloc(n * sizeof(float));
    if (argc == 2) {
        read_matr_formula(A, b, n);
    } else if (read_matr_file(A, b, n, argv[2]) != 0){
        printf("Can't read matrix from file");
        free(A);
        free(b);
        return 0;
    }

    float* x = (float*)malloc(n * sizeof(float));
    clock_t start = clock();
    if (solve_linear_system(A, n, b, x) != 0) {
        printf("Can't solve linear system");
        free(A);
        free(b);
        free(x);
        return 0;
    }
    clock_t end = clock();

    printf("Solved linear system in %f microseconds\n",
           ((double) (end - start)) / CLOCKS_PER_SEC * 1000);
    print_solution(x, n, MAX_PRINT_SIZE);
    print_error(A, n, b, x);
    if (argc == 2) {
        print_solution_error(x, n);
    }

    free(A);
    free(x);
    free(b);
    return 0;
}