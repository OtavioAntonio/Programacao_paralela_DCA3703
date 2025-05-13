#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#define NX 32
#define NY 32
#define NZ 32
#define NSTEPS 100000
#define DX 1.0
#define DY 1.0
#define DZ 1.0
#define DT 0.005
#define NU 0.05
#define SAVE_EVERY 20000

int main() {
    double u[NX][NY][NZ], v[NX][NY][NZ], w[NX][NY][NZ];
    double u_new[NX][NY][NZ], v_new[NX][NY][NZ], w_new[NX][NY][NZ];

    // Inicialização
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                u[i][j][k] = 0.0;
                v[i][j][k] = 0.0;
                w[i][j][k] = 0.0;
            }
        }
    }

    // Perturbação inicial no centro
    int cx = NX / 2, cy = NY / 2, cz = NZ / 2;
    u[cx][cy][cz] = 1.0;

    struct timeval start, end;
    gettimeofday(&start, NULL);  // Início da medição

    for (int n = 0; n <= NSTEPS; n++) {

        // Atualização com difusão (cálculo do laplaciano embutido)
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    double lap_u = (u[i+1][j][k] - 2*u[i][j][k] + u[i-1][j][k]) / (DX*DX)
                                 + (u[i][j+1][k] - 2*u[i][j][k] + u[i][j-1][k]) / (DY*DY)
                                 + (u[i][j][k+1] - 2*u[i][j][k] + u[i][j][k-1]) / (DZ*DZ);

                    double lap_v = (v[i+1][j][k] - 2*v[i][j][k] + v[i-1][j][k]) / (DX*DX)
                                 + (v[i][j+1][k] - 2*v[i][j][k] + v[i][j-1][k]) / (DY*DY)
                                 + (v[i][j][k+1] - 2*v[i][j][k] + v[i][j][k-1]) / (DZ*DZ);

                    double lap_w = (w[i+1][j][k] - 2*w[i][j][k] + w[i-1][j][k]) / (DX*DX)
                                 + (w[i][j+1][k] - 2*w[i][j][k] + w[i][j-1][k]) / (DY*DY)
                                 + (w[i][j][k+1] - 2*w[i][j][k] + w[i][j][k-1]) / (DZ*DZ);

                    u_new[i][j][k] = u[i][j][k] + NU * DT * lap_u;
                    v_new[i][j][k] = v[i][j][k] + NU * DT * lap_v;
                    w_new[i][j][k] = w[i][j][k] + NU * DT * lap_w;
                }
            }
        }

        // Cópia dos valores atualizados
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    u[i][j][k] = u_new[i][j][k];
                    v[i][j][k] = v_new[i][j][k];
                    w[i][j][k] = w_new[i][j][k];
                }
            }
        }

        gettimeofday(&end, NULL);  // Fim da medição

        // Salva os dados em CSV (opcional)
        
        if (n % SAVE_EVERY == 0) {
            char filename[128];
            snprintf(filename, sizeof(filename), "velocidade_u_fatia_z_t%04d.csv", n);
            FILE* fp = fopen(filename, "w");
            if (!fp) {
                perror("Erro ao abrir o arquivo");
                exit(1);
            }

            int z = NZ / 2;
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    fprintf(fp, "%.5f", u[i][j][z]);
                    if (j < NY - 1) fprintf(fp, ",");
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("Salvo: %s\n", filename);
        }
    
    }

    double tempo = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
    printf("Tempo de execução: %.6f segundos\n", tempo);

    return 0;
}


/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>  // Necessário para gettimeofday

#define NX 32
#define NY 32
#define NZ 32
#define NSTEPS 100000
#define DX 1.0
#define DY 1.0
#define DZ 1.0
#define DT 0.005
#define NU 0.05
#define SAVE_EVERY 20000  // Salvar a cada 100 passos

double u[NX][NY][NZ], v[NX][NY][NZ], w[NX][NY][NZ];
double u_new[NX][NY][NZ], v_new[NX][NY][NZ], w_new[NX][NY][NZ];

// Calcula o laplaciano
double laplacian(double field[NX][NY][NZ], int i, int j, int k) {
    double lap = 0.0;
    if (i > 0 && i < NX-1 && j > 0 && j < NY-1 && k > 0 && k < NZ-1) {
        lap = (field[i+1][j][k] - 2*field[i][j][k] + field[i-1][j][k]) / (DX*DX)
            + (field[i][j+1][k] - 2*field[i][j][k] + field[i][j-1][k]) / (DY*DY)
            + (field[i][j][k+1] - 2*field[i][j][k] + field[i][j][k-1]) / (DZ*DZ);
    }
    return lap;
}

void step() {
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
            for (int k = 1; k < NZ-1; k++) {
                u_new[i][j][k] = u[i][j][k] + NU * DT * laplacian(u, i, j, k);
                v_new[i][j][k] = v[i][j][k] + NU * DT * laplacian(v, i, j, k);
                w_new[i][j][k] = w[i][j][k] + NU * DT * laplacian(w, i, j, k);
            }
        }
    }

    // Atualiza as variáveis
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
            for (int k = 1; k < NZ-1; k++) {
                u[i][j][k] = u_new[i][j][k];
                v[i][j][k] = v_new[i][j][k];
                w[i][j][k] = w_new[i][j][k];
            }
        }
    }
}

void initialize() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                u[i][j][k] = 0.0;
                v[i][j][k] = 0.0;
                w[i][j][k] = 0.0;
            }
        }
    }

    // Perturbação inicial no centro
    int cx = NX / 2, cy = NY / 2, cz = NZ / 2;
    u[cx][cy][cz] = 1.0;
}

void save_to_csv_time(const char* prefix, int t) {
    char filename[128];
    snprintf(filename, sizeof(filename), "%s_t%04d.csv", prefix, t);

    FILE* fp = fopen(filename, "w");
    if (!fp) {
        perror("Erro ao abrir o arquivo CSV");
        exit(1);
    }

    int z = NZ / 2;
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            fprintf(fp, "%.5f", u[i][j][z]);
            if (j < NY - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Salvo: %s\n", filename);
}

int main() {
    // Definir o número de threads OpenMP
    //omp_set_num_threads(4); // Ajuste o número de threads conforme necessário

    struct timeval start, end;

    initialize();

    // Inicia a medição de tempo
    gettimeofday(&start, NULL);

    for (int n = 0; n <= NSTEPS; n++) {
        step();

        // Salvar fatias ao longo do tempo
        if (n % SAVE_EVERY == 0) {
            save_to_csv_time("velocidade_u_fatia_z", n);
        }
    }

    // Finaliza a medição de tempo
    gettimeofday(&end, NULL);

    // Calcula a diferença entre os tempos
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Tempo de execução: %.6f segundos\n", time_taken);

    return 0;
}
*/

/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>  // Necessário para gettimeofday

#define NX 32
#define NY 32
#define NZ 32
#define NSTEPS 1000
#define DX 1.0
#define DY 1.0
#define DZ 1.0
#define DT 0.005
#define NU 0.05

double u[NX][NY][NZ], v[NX][NY][NZ], w[NX][NY][NZ];
double u_new[NX][NY][NZ], v_new[NX][NY][NZ], w_new[NX][NY][NZ];

// Calcula o laplaciano
double laplacian(double field[NX][NY][NZ], int i, int j, int k) {
    double lap = 0.0;
    if (i > 0 && i < NX-1 && j > 0 && j < NY-1 && k > 0 && k < NZ-1) {
        lap = (field[i+1][j][k] - 2*field[i][j][k] + field[i-1][j][k]) / (DX*DX)
            + (field[i][j+1][k] - 2*field[i][j][k] + field[i][j-1][k]) / (DY*DY)
            + (field[i][j][k+1] - 2*field[i][j][k] + field[i][j][k-1]) / (DZ*DZ);
    }
    return lap;
}

void step() {
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
            for (int k = 1; k < NZ-1; k++) {
                u_new[i][j][k] = u[i][j][k] + NU * DT * laplacian(u, i, j, k);
                v_new[i][j][k] = v[i][j][k] + NU * DT * laplacian(v, i, j, k);
                w_new[i][j][k] = w[i][j][k] + NU * DT * laplacian(w, i, j, k);
            }
        }
    }

    // Atualiza as variáveis
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
            for (int k = 1; k < NZ-1; k++) {
                u[i][j][k] = u_new[i][j][k];
                v[i][j][k] = v_new[i][j][k];
                w[i][j][k] = w_new[i][j][k];
            }
        }
    }
}

void initialize() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                u[i][j][k] = 0.0;
                v[i][j][k] = 0.0;
                w[i][j][k] = 0.0;
            }
        }
    }

    // Perturbação inicial no centro
    int cx = NX / 2, cy = NY / 2, cz = NZ / 2;
    u[cx][cy][cz] = 1.0;
}

void save_to_csv(const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        perror("Erro ao abrir o arquivo CSV");
        exit(1);
    }

    int z = NZ / 2;
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            fprintf(fp, "%.5f", u[i][j][z]);
            if (j < NY - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Resultados salvos em %s\n", filename);
}

int main() {
    // Definir o número de threads OpenMP
    omp_set_num_threads(4); // Ajuste o número de threads conforme necessário

    struct timeval start, end;

    initialize();

    // Inicia a medição de tempo
    gettimeofday(&start, NULL);

    for (int n = 0; n < NSTEPS; n++) {
        step();
    }
    
    // Finaliza a medição de tempo
    gettimeofday(&end, NULL);
    
    save_to_csv("velocidade_u_fatia_z_paralelo.csv");

    

    // Calcula a diferença entre os tempos
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("Tempo de execução: %.6f segundos\n", time_taken);

    return 0;
}
*/