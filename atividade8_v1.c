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

    srand(time(NULL)); // Semente global (não thread-safe)

    struct timeval start, end;

    gettimeofday(&start, NULL);

    #pragma omp parallel private(privateHit)
    {
        privateHit = 0;

        #pragma omp for
        for (long long i = 0; i <= n; i++) {
            double x, y;

            
            x = ((double)rand() / RAND_MAX) * l;
            y = ((double)rand() / RAND_MAX) * l;
           

            if (x * x + y * y <= 1) {
                privateHit++;
            }
        }

        #pragma omp critical
        {
            hit += privateHit;
        }
    }

    gettimeofday(&end, NULL);
    double tempo_execucao = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;

    pi = (l * l * hit) / (r * r * n);
    printf("Tempo de execução: %.6f segundos\n", tempo_execucao);
    printf("Pi: %.6f\n", pi);

    return 0;
}
