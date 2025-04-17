#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>

//#define SIZE 100000000 // 100 milhões de elementos

int main() {

    long long n, cont = 0, primo = 0;

    printf("Digite n: ");
    scanf("%lld", &n);

    // Variáveis para medir o tempo
    struct timeval start, end;
    
    // Obtém o tempo de início
    gettimeofday(&start, NULL);

    for (int i = 2; i <= n; i++) {
        cont = 0;
        //printf("i: %d\n", i);
        for(int j=1; j<=i; j++){
            //printf("j: %d\n", j);
            if(i%j==0) cont++;
        }
       
        if(cont==2) primo++;
        
    }

    // Obtém o tempo de fim
    gettimeofday(&end, NULL);
    
    // Calcula o tempo de execução em segundos e microssegundos
    double tempo_execucao = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    
    // Exibe o resultado
    printf("Tempo de execução Sequencial: %.6f segundos\n", tempo_execucao);
    printf("Qtde primos: %lld\n", primo);

    primo = 0;

    // Obtém o tempo de início
    gettimeofday(&start, NULL);

    //#pragma omp parallel
    #pragma omp parallel for private(cont) reduction(+:primo)
        for (int i = 2; i <= n; i++) {
            cont = 0;
            //printf("i: %d\n", i);
            for(int j=1; j<=i; j++){
                //printf("j: %d\n", j);
                if(i%j==0) cont++;
            }
        
            if(cont==2) primo++;
            
        }
    // Obtém o tempo de fim
    gettimeofday(&end, NULL);
    
    // Calcula o tempo de execução em segundos e microssegundos
    tempo_execucao = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    
    // Exibe o resultado
    printf("Tempo de execução Paralelo: %.6f segundos\n", tempo_execucao);
    printf("Qtde primos: %lld\n", primo);

    return 0;
}
