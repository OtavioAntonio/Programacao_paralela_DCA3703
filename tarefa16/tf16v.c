#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void matrix_vector_multiply(double *A, double *x, double *y, int start_row, int end_row, int N) {
    for (int i = start_row; i < end_row; i++) {
        y[i] = 0.0;
        for (int j = 0; j < N; j++) {
            y[i] += A[i * N + j] * x[j];
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int start_size = 50;  // tamanho inicial da matriz
    int max_size = 50;    // tamanho máximo
    int step = 10;        // passo do loop

    for (int M = start_size; M <= max_size; M += step) {
        int N = M;  // matriz quadrada

        if (M % size != 0) {
            if (rank == 0) {
                printf("Erro: M (%d) deve ser divisível pelo número de processos (%d).\n", M, size);
            }
            MPI_Finalize();
            return -1;
        }

        double *A = NULL;
        double *x = (double *)malloc(N * sizeof(double));
        double *y = (double *)malloc(M * sizeof(double));

        if (rank == 0) {
            A = (double *)malloc(M * N * sizeof(double));
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    A[i * N + j] = i + j;
                }
            }

            for (int i = 0; i < N; i++) {
                x[i] = 1.0;
            }
        }

        

        int linhas_per_proc = M / size;

        double *local_A = (double *)malloc(linhas_per_proc * N * sizeof(double));
        double *local_y = (double *)malloc(linhas_per_proc * sizeof(double));

        MPI_Barrier(MPI_COMM_WORLD); // sincroniza todos os processos

        double start = MPI_Wtime();

        MPI_Scatter(A, linhas_per_proc * N, MPI_DOUBLE,
                    local_A, linhas_per_proc * N, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        matrix_vector_multiply(local_A, x, local_y, 0, linhas_per_proc, N);

        MPI_Gather(local_y, linhas_per_proc, MPI_DOUBLE,
                   y, linhas_per_proc, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);

        double end = MPI_Wtime();

        if (rank == 0) {
            printf("Tempo de execução para matriz %dx%d: %lf segundos\n", M, N, end - start);
        }

        free(A);
        free(x);
        free(y);
        free(local_A);
        free(local_y);
    }

    MPI_Finalize();
    return 0;
}