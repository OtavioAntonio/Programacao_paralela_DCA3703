#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int rank; // id do processo atual
    int total_de_processos; // nº total de processos

    MPI_Init(&argc, &argv); // Inicializa o ambiente MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Obtém o id (rank) do processo atual
    MPI_Comm_size(MPI_COMM_WORLD, &total_de_processos); // Obtém o número total de processos

    // Verifica se o tamanho da matriz foi passado como argumento
    if (argc < 2) {
        if (rank == 0) {
            printf("Uso: %s <tamanho da matriz>\n", argv[0]);
        }
        MPI_Finalize();
        return -1;
    }
    // **********************************************************
    
    int qtd_linhas = atoi(argv[1]);
    int qtd_colunas = qtd_linhas;  // matriz quadrada

    // Verifica o tamanho da Matriz é divisível pelo número de processos
    if (qtd_colunas % total_de_processos != 0) {
        if (rank == 0) {
            printf("Erro: O tamanho da matriz (%d) deve ser divisível pelo número de processos (%d).\n", N, size);
        }
        MPI_Finalize();
        return -1;
    }
    // ************************************************

    int colunas_por_processo = qtd_colunas / total_de_processos;

    double *matriz_A = NULL; // matriz completa (apenas no rank 0)
    double *vetor_X = NULL; // vetor completo (apenas no rank 0)
    double *vetor_R = (double *)calloc(qtd_colunas, sizeof(double)); // vetor resultado final (apenas no rank 0)

    // Vetores locais para cada processo
    double *bloco_local_de_A = (double *)malloc(qtd_colunas * colunas_por_processo * sizeof(double));
    // bloco local de colunas
    double *bloco_local_de_X = (double *)malloc(colunas_por_processo * sizeof(double));
    // parte local do vetor x
    double *bloco_local_de_R = (double *)calloc(qtd_colunas, sizeof(double));
    // resultado local

    // Inicialização apenas no processo 0
    if (rank == 0) {
        matriz_A = (double *)malloc(qtd_colunas * qtd_linhas * sizeof(double)); // matriz completa
        vetor_X = (double *)malloc(qtd_linhas * sizeof(double)); // vetor completo

        // Inicializa a matriz_A (exemplo: A[i][j] = i + j). 
        // É armazenada assim pois não sabemos quantos elementos serão inseridos
        for (int i = 0; i < qtd_linhas; i++) {
            for (int j = 0; j < qtd_colunas; j++) {
                matriz_A[i * qtd_colunas + j] = i + j;
            }
        }

        // Inicializa o vetor x com 1.0
        for (int i = 0; i < qtd_colunas; i++) {
            x[i] = 1.0;
        }
    }

     // Criação de tipo derivado para representar blocos de colunas
    MPI_Datatype tipo_coluna, tipo_coluna_redimensionado;
    MPI_Type_vector(qtd_linhas, colunas_por_processo, qtd_colunas, MPI_DOUBLE, &tipo_coluna);
    // cria um tipo derivado que representa varios blocos de dados
    // qtd de linhas => quantidade de blocos
    // colunas_por_processo => quantidade de elementos em cada bloco
    // qtd_colunas => Distância entre o inicio de dois blocos consecutivos
    // MPI_double => tipo do dado
    // tipo_coluna => onde será gravado o tipo derivado criado
    MPI_Type_create_resized(tipo_coluna, 0, colunas_por_processo * sizeof(double), &tipo_coluna_redimensionado);
    // Ele só diz ao MPI que, depois de enviar/receber um desses
    // tipos derivados, o próximo deveria começar 8 bytes depois 
    MPI_Type_commit(&tipo_coluna_redimensionado);

    // Distribui as colunas da matriz A entre os processos
    MPI_Scatter(matriz_A, 1, tipo_coluna_redimensionado, bloco_local_de_A, qtd_linhas * colunas_por_processo, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    /// ****** abaixo não está editado ***********


    // Distribui as partes do vetor x correspondentes às colunas
    if (rank == 0) {
        MPI_Scatter(x, cols_per_proc, MPI_DOUBLE,
                    local_x, cols_per_proc, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
    } else {
        MPI_Scatter(NULL, cols_per_proc, MPI_DOUBLE,
                    local_x, cols_per_proc, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
    }

    // Sincroniza antes de medir o tempo
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // Cálculo da contribuição parcial: local_y = A_local * x_local
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < cols_per_proc; j++) {
            local_y[i] += local_A[i * cols_per_proc + j] * local_x[j];
        }
    }

    // Reduz as contribuições locais no vetor y final no rank 0
    MPI_Reduce(local_y, y, M, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double end = MPI_Wtime();

    if (rank == 0) {
        printf("Tempo de execução para matriz %dx%d: %lf segundos\n", M, N, end - start);
        // Exibe parte do resultado (opcional)
        printf("Primeiros 5 elementos do vetor y:\n");
        for (int i = 0; i < (M < 5 ? M : 5); i++) {
            printf("y[%d] = %lf\n", i, y[i]);
        }
    }

    // Libera memória
    if (rank == 0) {
        free(A);
        free(x);
        free(y);
    }

    free(local_A);
    free(local_x);
    free(local_y);
    MPI_Type_free(&column_type);
    MPI_Type_free(&resized_col_type);

    MPI_Finalize();
    return 0;
}