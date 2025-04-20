#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <math.h>

int main() {
    long long n, hit = 0, privateHit = 0;
    double l = 1.0, r, pi;
    r = l / 2.0;

    printf("Digite n: ");
    scanf("%lld", &n);

    // Primeiro cálculo: hit compartilhado (sem private)
    hit = 0;

    // Variáveis para medir o tempo
    struct timeval start, end;
    
    // Obtém o tempo de início
    gettimeofday(&start, NULL);


    #pragma omp parallel for
    for (long long i = 0; i <= n; i++) {
        unsigned int seed = 1234 + omp_get_thread_num() * 1000 + i;
        double x = ((double)rand_r(&seed) / RAND_MAX) * l;
        double y = ((double)rand_r(&seed) / RAND_MAX) * l;

        if (x * x + y * y <= 1) {
            #pragma omp critical
            hit++;
        }
    }

    // Obtém o tempo de fim
    gettimeofday(&end, NULL);
    
    // Calcula o tempo de execução em segundos e microssegundos
    double tempo_execucao = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    // Exibe o resultado
    printf("Tempo de execução Sequencial: %.6f segundos\n", tempo_execucao);
    pi = (l * l * hit) / (r * r * n);
    printf("Pi (variáveis compartilhadas): %.6f\n", pi);



    // Segundo cálculo: hit privado por thread + critical
    pi = 0;
    hit = 0;

    // Obtém o tempo de início
    gettimeofday(&start, NULL);

    #pragma omp parallel private(privateHit)
    {
        privateHit = 0;

        #pragma omp for
        for (long long i = 0; i <= n; i++) {
            unsigned int seed = 1234 + omp_get_thread_num() * 1000 + i;
            double x = ((double)rand_r(&seed) / RAND_MAX) * l;
            double y = ((double)rand_r(&seed) / RAND_MAX) * l;

            if (x * x + y * y <= 1) {
                privateHit++;
            }
        }

        #pragma omp critical
        {
            hit += privateHit;
        }
    }

    // Obtém o tempo de fim
    gettimeofday(&end, NULL);
    
    // Calcula o tempo de execução em segundos e microssegundos
    tempo_execucao = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    // Exibe o resultado
    printf("Tempo de execução Sequencial: %.6f segundos\n", tempo_execucao);
    pi = (l * l * hit) / (r * r * n);
    printf("Pi (variáveis privadas): %.6f\n", pi);

    return 0;
}
