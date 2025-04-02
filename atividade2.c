#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N 100000000 // Tamanho do vetor

double get_tempo(struct timeval start, struct timeval end) {
    return (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
}

int main() {
    int *vetor = (int*)malloc(N * sizeof(int));
    int a = 2; // Fator de multiplicação
    struct timeval start, end;
    double tempo_decorrido;
    long long soma = 0, soma_otimizada = 0;
    long long temp1 = 0, temp2 = 0;

    // Preenchendo o vetor
    gettimeofday(&start, NULL);
    for (int i = 0; i < N; i++) {
        vetor[i] = a + i;
    }
    gettimeofday(&end, NULL);
    printf("Tempo de execução para preencher o vetor: %f segundos\n", get_tempo(start, end));

    // Soma acumulada
    gettimeofday(&start, NULL);
    for (int i = 0; i < N; i++) {
        soma += vetor[i];
    }
    gettimeofday(&end, NULL);
    printf("Tempo de execução para soma acumulada: %f segundos\n", get_tempo(start, end));
    printf("Resultado soma acumulada: %lld\n", soma);

    // Soma acumulada otimizada
    gettimeofday(&start, NULL);
    for (int i = 0; i < N; i+=2) {
        temp1 += vetor[i];
        temp2 += vetor[i+1];
    }
    soma_otimizada = temp1 + temp2;
    gettimeofday(&end, NULL);
    printf("Tempo de execução para soma acumulada otimizada: %f segundos\n", get_tempo(start, end));
    printf("Resultado soma otimizada: %lld\n", soma_otimizada);

    free(vetor);
    return 0;
}