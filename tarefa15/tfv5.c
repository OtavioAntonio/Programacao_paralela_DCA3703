#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 1000000     // Total de pontos na barra
#define STEPS 1000     // Número de iterações
#define LEFT 0
#define RIGHT 1

void initialize(double *local, int local_n, int rank) {
    for (int i = 1; i <= local_n; i++) {
        local[i] = rank;  // valor inicial depende do rank
    }
}

void simulate_blocking(double *local, double *temp, int local_n, int rank, int size) {
    int left = rank - 1;
    int right = rank + 1;

    for (int step = 0; step < STEPS; step++) {
        // Troca de bordas
        if (left >= 0)
            MPI_Send(&local[1], 1, MPI_DOUBLE, left, RIGHT, MPI_COMM_WORLD);
        if (right < size)
            MPI_Send(&local[local_n], 1, MPI_DOUBLE, right, LEFT, MPI_COMM_WORLD);

        if (left >= 0)
            MPI_Recv(&local[0], 1, MPI_DOUBLE, left, LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (right < size)
            MPI_Recv(&local[local_n + 1], 1, MPI_DOUBLE, right, RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Atualiza os pontos internos
        for (int i = 1; i <= local_n; i++) {
            temp[i] = 0.25 * (local[i - 1] + 2 * local[i] + local[i + 1]);
        }

        // Copia temp para local
        memcpy(&local[1], &temp[1], local_n * sizeof(double));
    }
}

void simulate_nonblocking(double *local, double *temp, int local_n, int rank, int size) {
    int left = rank - 1;
    int right = rank + 1;

    MPI_Request reqs[4]; // [0] recv_esq, [1] recv_dir, [2] send_esq, [3] send_dir

    for (int step = 0; step < STEPS; step++) {
        // Comunicação não bloqueante: Irecv primeiro
        if (left >= 0) MPI_Irecv(&local[0], 1, MPI_DOUBLE, left, LEFT, MPI_COMM_WORLD, &reqs[0]);
        else reqs[0] = MPI_REQUEST_NULL;

        if (right < size) MPI_Irecv(&local[local_n + 1], 1, MPI_DOUBLE, right, RIGHT, MPI_COMM_WORLD, &reqs[1]);
        else reqs[1] = MPI_REQUEST_NULL;

        // Isend depois
        if (left >= 0) MPI_Isend(&local[1], 1, MPI_DOUBLE, left, RIGHT, MPI_COMM_WORLD, &reqs[2]);
        else reqs[2] = MPI_REQUEST_NULL;

        if (right < size) MPI_Isend(&local[local_n], 1, MPI_DOUBLE, right, LEFT, MPI_COMM_WORLD, &reqs[3]);
        else reqs[3] = MPI_REQUEST_NULL;

        // Atualiza os pontos internos (exceto bordas)
        for (int i = 2; i <= local_n - 1; i++) {
            temp[i] = 0.25 * (local[i - 1] + 2 * local[i] + local[i + 1]);
        }

        // Espera recepção da esquerda, se existir, e atualiza a borda
        if (left >= 0) {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
            temp[1] = 0.25 * (local[0] + 2 * local[1] + local[2]);
        } else temp[1] = 0.25 * (local[1] + 2 * local[1] + local[2]); // condição de contorno fixa

        // Espera recepção da direita, se existir, e atualiza a borda
        if (right < size) {
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
            temp[local_n] = 0.25 * (local[local_n - 1] + 2 * local[local_n] + local[local_n + 1]);
        } else temp[local_n] = 0.25 * (local[local_n - 1] + 2 * local[local_n] + local[local_n]); // condição de contorno fixa

        // Espera os envios finalizarem (boa prática)
        if (left >= 0) MPI_Wait(&reqs[2], MPI_STATUS_IGNORE);
        if (right < size) MPI_Wait(&reqs[3], MPI_STATUS_IGNORE);

        // Copia o estado atualizado para o vetor principal
        memcpy(&local[1], &temp[1], local_n * sizeof(double));
    }
}

void simulate_overlap(double *local, double *temp, int local_n, int rank, int size) {
    int left = rank - 1;
    int right = rank + 1;

    MPI_Request reqs[4];
    int flag[2] = {0, 0};

    for (int step = 0; step < STEPS; step++) {
        // Inicia comunicações
        if (left >= 0)
            MPI_Irecv(&local[0], 1, MPI_DOUBLE, left, LEFT, MPI_COMM_WORLD, &reqs[0]);
        if (right < size)
            MPI_Irecv(&local[local_n + 1], 1, MPI_DOUBLE, right, RIGHT, MPI_COMM_WORLD, &reqs[1]);
        if (left >= 0)
            MPI_Isend(&local[1], 1, MPI_DOUBLE, left, RIGHT, MPI_COMM_WORLD, &reqs[2]);
        if (right < size)
            MPI_Isend(&local[local_n], 1, MPI_DOUBLE, right, LEFT, MPI_COMM_WORLD, &reqs[3]);

        // Computa os pontos internos que não dependem da borda
        for (int i = 2; i <= local_n - 1; i++) {
            temp[i] = 0.25 * (local[i - 1] + 2 * local[i] + local[i + 1]);
        }

        // Aguarda comunicação da esquerda e atualiza
        if (left >= 0) {
            while (!flag[0])
                MPI_Test(&reqs[0], &flag[0], MPI_STATUS_IGNORE);
            temp[1] = 0.25 * (local[0] + 2 * local[1] + local[2]);
        } else {
            temp[1] = 0.25 * (local[1] + 2 * local[1] + local[2]);  // fronteira fixa
        }

        // Aguarda comunicação da direita e atualiza
        if (right < size) {
            while (!flag[1])
                MPI_Test(&reqs[1], &flag[1], MPI_STATUS_IGNORE);
            temp[local_n] = 0.25 * (local[local_n - 1] + 2 * local[local_n] + local[local_n + 1]);
        } else {
            temp[local_n] = 0.25 * (local[local_n - 1] + 2 * local[local_n] + local[local_n]);  // fronteira fixa
        }

        // Espera envios completarem
        if (left >= 0) MPI_Wait(&reqs[2], MPI_STATUS_IGNORE);
        if (right < size) MPI_Wait(&reqs[3], MPI_STATUS_IGNORE);

        // Atualiza estado
        memcpy(&local[1], &temp[1], local_n * sizeof(double));

        // Reset flags
        flag[0] = flag[1] = 0;
    }
}

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_n = N / size;
    double *local = malloc((local_n + 2) * sizeof(double));  // inclui células fantasmas
    double *temp = malloc((local_n + 2) * sizeof(double));

    initialize(local, local_n, rank);

    double start, end;

    if (rank == 0) printf("Iniciando simulação com comunicação bloqueante...\n");
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    simulate_blocking(local, temp, local_n, rank, size);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if (rank == 0) printf("Tempo (bloqueante): %f segundos\n", end - start);

    initialize(local, local_n, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) printf("Iniciando simulação com comunicação não bloqueante...\n");
    start = MPI_Wtime();
    simulate_nonblocking(local, temp, local_n, rank, size);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if (rank == 0) printf("Tempo (não bloqueante): %f segundos\n", end - start);

    initialize(local, local_n, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) printf("Iniciando simulação com sobreposição comunicação-computação...\n");
    start = MPI_Wtime();
    simulate_overlap(local, temp, local_n, rank, size);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if (rank == 0) printf("Tempo (overlap): %f segundos\n", end - start);

    free(local);
    free(temp);
    MPI_Finalize();
    return 0;
}