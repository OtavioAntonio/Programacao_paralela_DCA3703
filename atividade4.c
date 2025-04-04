#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>

#define SIZE 100000000 // 100 milhões de elementos

int main() {
    int *a = (int *)malloc(SIZE * sizeof(int));
    double *b = (double *)malloc(SIZE * sizeof(double));

    if (a == NULL || b == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    int threads[] = {1, 2, 4, 8, 16, 32}; // Variação do número de threads

    for (int t = 0; t < 6; t++) {
        int num_threads = threads[t];
        printf("Starting parallel tasks with %d threads...\n", num_threads);

        struct timeval start, end;

        // Tarefa limitada por memória
        omp_set_num_threads(num_threads);
        gettimeofday(&start, NULL);
        #pragma omp parallel for
        for (int i = 0; i < SIZE; i++) {
            a[i] = i + 5;
        }
        gettimeofday(&end, NULL);
        double elapsed = (end.tv_sec - start.tv_sec) + 
                         (end.tv_usec - start.tv_usec) / 1e6;
        printf("Memory-bound task with %d threads took %.6f seconds.\n", num_threads, elapsed);

        // Tarefa limitada por CPU
        gettimeofday(&start, NULL);
        #pragma omp parallel for
        for (int i = 0; i < SIZE; i++) {
            b[i] = sqrt((i * 0.001) * cos(i * 0.001) + log(i + 1));
        }
        gettimeofday(&end, NULL);
        elapsed = (end.tv_sec - start.tv_sec) + 
                  (end.tv_usec - start.tv_usec) / 1e6;
        printf("CPU-bound task with %d threads took %.6f seconds.\n", num_threads, elapsed);

        printf("-----------------------------------------\n");
    }

    free(a);
    free(b);
}
