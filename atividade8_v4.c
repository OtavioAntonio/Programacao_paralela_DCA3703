#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <math.h>

int main() {
    long long n;
    double l = 1.0, r, pi;

    r = l / 2.0;

    printf("Digite n: ");
    scanf("%lld", &n);

    srand(time(NULL)); // Semente global (não thread-safe)

    struct timeval start, end;

    // Vetor para armazenar o número de acertos de cada thread
    long long *hits_por_thread = (long long *)malloc(sizeof(long long) * omp_get_max_threads());

    gettimeofday(&start, NULL);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num(); // Identificador da thread
        hits_por_thread[tid] = 0; // Inicializa o contador privado da thread
        unsigned int seed = time(NULL) + tid; // Cada thread cria sua própria semente

        #pragma omp for
        for (long long i = 0; i < n; i++) {
            double x, y;

            // Usa rand_r() que é thread-safe
            x = ((double)rand_r(&seed) / RAND_MAX) * l;
            y = ((double)rand_r(&seed) / RAND_MAX) * l;

            // Verifica se o ponto está dentro do círculo
            if (x * x + y * y <= 1) {
                hits_por_thread[tid]++; // Incrementa o contador da thread
            }
        }
    }
    // Obtém o tempo de fim
    gettimeofday(&end, NULL);

    // Após a execução paralela, soma os acertos de todas as threads
    long long hit_total = 0;
    for (int i = 0; i < omp_get_max_threads(); i++) {
        hit_total += hits_por_thread[i];
    }
    
    double tempo_execucao = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;

    // Calcula a estimativa de pi
    pi = (l * l * hit_total) / (r * r * n);

    printf("Tempo de execução: %.6f segundos\n", tempo_execucao);
    printf("Pi: %.6f\n", pi);

    // Libera a memória alocada para o vetor de hits
    free(hits_por_thread);

    return 0;
}