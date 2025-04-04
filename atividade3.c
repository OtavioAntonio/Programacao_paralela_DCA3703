#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>  // Para usar M_PI

int main() {
    int iteracoes;
    
    // Lê o número de iterações do usuário
    printf("Digite o número de iterações para calcular pi: ");
    scanf("%d", &iteracoes);
    
    // Variáveis para medir o tempo
    struct timeval start, end;
    
    // Obtém o tempo de início
    gettimeofday(&start, NULL);
    
    // Calcula a aproximação de pi usando a série de Leibniz
    double pi = 0.0;
    for (int i = 0; i < iteracoes; i++) {
        if (i % 2 == 0) {
            pi += 4.0 / (2 * i + 1);  // Termos positivos
        } else {
            pi -= 4.0 / (2 * i + 1);  // Termos negativos
        }
    }
    
    // Obtém o tempo de fim
    gettimeofday(&end, NULL);
    
    // Calcula o tempo de execução em segundos e microssegundos
    double tempo_execucao = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    
    // Exibe o resultado
    printf("Aproximação de Pi após %d iterações: %.20f\n", iteracoes, pi);
    printf("Valor de Pi da biblioteca (M_PI): %.20f\n", M_PI);
    printf("Tempo de execução: %.6f segundos\n", tempo_execucao);
    
    // Comparação com o valor de Pi da biblioteca
    double erro = fabs(pi - M_PI); // Erro absoluto
    printf("Erro absoluto entre a aproximação e o valor de Pi: %.20f\n", erro);
    
    return 0;
}