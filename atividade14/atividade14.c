#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_EXCHANGES 10000  // número de trocas para cada tamanho

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2) {
        if (rank == 0)
            printf("Este programa requer exatamente 2 processos.\n");
        MPI_Finalize();
        return 1;
    }

    // Tamanhos de mensagem em bytes (de 8 B a 1 MB)
    int sizes[] = {8, 64, 512, 4096, 32768, 262144, 1048576}; // até 1MB
    int num_sizes = sizeof(sizes) / sizeof(int);

    for (int i = 0; i < num_sizes; i++) {
        int msg_size = sizes[i];
        char *buffer = (char *)malloc(msg_size);
        memset(buffer, 'a', msg_size);  // preenche o buffer com dados

        MPI_Barrier(MPI_COMM_WORLD);  // sincroniza antes de medir o tempo
        double start_time = MPI_Wtime();

        for (int j = 0; j < NUM_EXCHANGES; j++) {
            if (rank == 0) {
                MPI_Send(buffer, msg_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(buffer, msg_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else if (rank == 1) {
                MPI_Recv(buffer, msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(buffer, msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }

        double end_time = MPI_Wtime();
        if (rank == 0) {
            double total_time = end_time - start_time;
            printf("Tamanho: %7d bytes | Tempo total: %8.6f s | Tempo médio por troca: %10.8f s\n",
                   msg_size, total_time, total_time / NUM_EXCHANGES);
        }

        free(buffer);
    }

    MPI_Finalize();
    return 0;
}
