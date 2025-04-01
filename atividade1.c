#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

void multiplicarMatrizVetorLECI(int **matriz, int *vetor, int *resultado, int linhas, int colunas) {
    for (int i = 0; i < linhas; i++) {
        resultado[i] = 0; // Inicializa o resultado da linha
        for (int j = 0; j < colunas; j++) {
            resultado[i] += matriz[i][j] * vetor[j]; // Multiplicação e soma
        }
    }
}
void multiplicarMatrizVetorCELI(int **matriz, int *vetor, int *resultado, int linhas, int colunas) {
    for (int j = 0; j < colunas; j++) {
        resultado[j] = 0; // Inicializa o resultado da linha
        for (int i = 0; i < linhas; i++) {
            resultado[j] += matriz[i][j] * vetor[i]; // Multiplicação e soma
        }
    }
}
int main() {
    int linhas, colunas;

    // Solicita ao usuário o número de linhas e colunas da matriz
    printf("Digite o n\u00famero de linhas da matriz: ");
    scanf("%d", &linhas);
    printf("Digite o n\u00famero de colunas da matriz: ");
    scanf("%d", &colunas);

    // Aloca dinamicamente a matriz
    int **matriz = (int **)malloc(linhas * sizeof(int *));
    for (int i = 0; i < linhas; i++) {
        matriz[i] = (int *)malloc(colunas * sizeof(int));
    }

    // Aloca dinamicamente o vetor
    int *vetor = (int *)malloc(colunas * sizeof(int));
    
    // Aloca dinamicamente o vetor resultado
    int *resultado = (int *)malloc(linhas * sizeof(int));

    // Garante aleatoriedade a cada execução
    srand(time(NULL));

    // Preenche a matriz com valores aleatórios
    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            matriz[i][j] = rand() % 100;
        }
    }

    // Preenche o vetor com valores aleatórios
    for (int i = 0; i < colunas; i++) {
        vetor[i] = rand() % 100;
    }

    // Marca o tempo antes de iniciar a execução do cálculo
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // Multiplica a matriz pelo vetor
    multiplicarMatrizVetorCELI(matriz, vetor, resultado, linhas, colunas);

    // Marca o tempo após a execução do cálculo
    gettimeofday(&end, NULL);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Tempo de execu\u00e7\u00e3o Linha interna - Coluna externa (Walltime): %.8f segundos\n", time_taken);

    // Marca o tempo antes de iniciar a execução do cálculo
    gettimeofday(&start, NULL);

    // Multiplica a matriz pelo vetor
    multiplicarMatrizVetorLECI(matriz, vetor, resultado, linhas, colunas);

    // Marca o tempo após a execução do cálculo
    gettimeofday(&end, NULL);

    time_taken = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Tempo de execu\u00e7\u00e3o Linha externa - Coluna interna (CONTIGUA) (Walltime): %.8f segundos\n", time_taken);

    // Libera a memória alocada
    for (int i = 0; i < linhas; i++) {
        free(matriz[i]);
    }
    free(matriz);
    free(vetor);
    free(resultado);

    return 0;
}