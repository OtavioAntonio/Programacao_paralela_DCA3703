#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int start_size = 10;  // tamanho inicial da matriz
    int max_size = 10;    // para teste mais rápido
    int step = 10;

    for (int M = start_size; M <= max_size; M += step) {
        int N = M;  // matriz quadrada

        if (N % size != 0) {
            if (rank == 0) {
                printf("Erro: N (%d) deve ser divisível pelo número de processos (%d).\n", N, size);
            }
            MPI_Finalize();
            return -1;
        }

        int cols_per_proc = N / size;

        double *A = NULL;
        double *x = NULL;
        double *y = (double *)calloc(M, sizeof(double));  // resultado final (no rank 0)
        double *local_A = (double *)malloc(M * cols_per_proc * sizeof(double)); // bloco local de colunas
        double *local_x = (double *)malloc(cols_per_proc * sizeof(double));     // parte de x local
        double *local_y = (double *)calloc(M, sizeof(double));                  // contribuição local

        if (rank == 0) {
            A = (double *)malloc(M * N * sizeof(double));
            x = (double *)malloc(N * sizeof(double));

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    A[i * N + j] = i + j; // Exemplo: A[i][j] = i + j
                }
            }

            for (int i = 0; i < N; i++) {
                x[i] = 1.0;
            }
        }

        // Criar tipo derivado para uma coluna
        MPI_Datatype column_type, resized_col_type;
        MPI_Type_vector(M, cols_per_proc, N, MPI_DOUBLE, &column_type);
        MPI_Type_create_resized(column_type, 0, cols_per_proc * sizeof(double), &resized_col_type);
        MPI_Type_commit(&resized_col_type);

        // Scatter das colunas da matriz A
        MPI_Scatter(A, 1, resized_col_type,
                    local_A, M * cols_per_proc, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        // Scatter dos pedaços do vetor x correspondentes às colunas
        if (rank == 0) {
            MPI_Scatter(x, cols_per_proc, MPI_DOUBLE,
                        local_x, cols_per_proc, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
        } else {
            MPI_Scatter(NULL, cols_per_proc, MPI_DOUBLE,
                        local_x, cols_per_proc, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
        }

        // Sincroniza antes da medição de tempo
        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        // Cálculo da contribuição parcial: y += A_local * x_local
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < cols_per_proc; j++) {
                local_y[i] += local_A[i * cols_per_proc + j] * local_x[j];
            }
        }

        // Reduz os vetores locais em y global no rank 0
        MPI_Reduce(local_y, y, M, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        double end = MPI_Wtime();

        if (rank == 0) {
            printf("Tempo de execução para matriz %dx%d: %lf segundos\n", M, N, end - start);
        }

        // Liberação
        if (rank == 0) {
            free(A);
            free(x);
            free(y);
        }

        free(local_A);
        free(local_x);
        free(local_y);
        MPI_Type_free(&column_type);
        MPI_Type_free(&resized_col_type);
    }

    MPI_Finalize();
    return 0;
}