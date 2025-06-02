#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int M = size;  // número de processos deve ser igual ao número de linhas
    int N = M;     // matriz quadrada MxN

    double *A = NULL;
    double *x = (double *)malloc(N * sizeof(double));
    double *y = NULL;

    if (rank == 0) {
        A = (double *)malloc(M * N * sizeof(double));
        y = (double *)malloc(M * sizeof(double));

        // inicializa matriz A com A[i][j] = i + j
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                A[i * N + j] = i + j;
            }
        }

        // inicializa vetor x com 1.0
        for (int i = 0; i < N; i++) {
            x[i] = 1.0;
        }
    }

    // broadcast do vetor x para todos os processos
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // cada processo recebe uma linha da matriz A
    double *local_A = (double *)malloc(N * sizeof(double));
    MPI_Scatter(A, N, MPI_DOUBLE,
                local_A, N, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // cada processo multiplica sua linha pelo vetor x
    double local_y = 0.0;
    for (int j = 0; j < N; j++) {
        local_y += local_A[j] * x[j];
    }

    // cada processo envia seu resultado (y[i]) para o processo 0
    MPI_Gather(&local_y, 1, MPI_DOUBLE,
               y, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Resultado y = A * x:\n");
        for (int i = 0; i < M; i++) {
            printf("y[%d] = %f\n", i, y[i]);
        }
    }

    free(local_A);
    free(x);
    if (rank == 0) {
        free(A);
        free(y);
    }

    MPI_Finalize();
    return 0;
}
