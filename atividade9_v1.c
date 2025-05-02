
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <unistd.h> // sleep()

#define N 10  // Poucas inserções para visualização

typedef struct Node {
    int value;
    struct Node* next;
} Node;

void insert(Node** head, int value) {
    Node* new_node = malloc(sizeof(Node));
    new_node->value = value;
    new_node->next = *head;
    *head = new_node;
}

void print_list(Node* head, const char* name) {
    printf("%s: ", name);
    Node* current = head;
    while (current) {
        printf("%d -> ", current->value);
        current = current->next;
    }
    printf("NULL\n");
}

void free_list(Node* head) {
    Node* tmp;
    while (head) {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

int main() {
    Node* lista1 = NULL;
    Node* lista2 = NULL;

    srand(time(NULL));

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 0; i < N; i++) {
                int value = rand() % 100;
                #pragma omp task firstprivate(value, i)
                {
                    int tid = omp_get_thread_num();
                    int list_choice = rand() % 2;

                    if (list_choice == 0) {
                        #pragma omp critical(lista1)  // experimente tirar o (lista1)
                        {
                            printf("Thread %d (iteracao %d): INSERINDO %d na Lista 1...\n", tid, i, value);
                            insert(&lista1, value);
                            printf("Thread %d (iteracao %d): FINALIZOU insercao na Lista 1\n", tid, i);
                        }
                    } else {
                        #pragma omp critical(lista2)  // experimente tirar o (lista2)
                        {
                            printf("Thread %d (iteracao %d): INSERINDO %d na Lista 2...\n", tid, i, value);
                            insert(&lista2, value);
                            printf("Thread %d (iteracao %d): FINALIZOU insercao na Lista 2\n", tid, i);
                        }
                    }
                }
            }
        }
    }

    printf("\nListas finais:\n");
    print_list(lista1, "Lista 1");
    print_list(lista2, "Lista 2");

    free_list(lista1);
    free_list(lista2);

    return 0;
}

/*
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <unistd.h> // sleep()

#define N 10

typedef struct Node {
    int value;
    struct Node* next;
} Node;

void insert(Node** head, int value) {
    Node* new_node = malloc(sizeof(Node));
    new_node->value = value;
    new_node->next = *head;
    *head = new_node;
}

void print_list(Node* head, const char* name) {
    printf("%s: ", name);
    Node* current = head;
    while (current) {
        printf("%d -> ", current->value);
        current = current->next;
    }
    printf("NULL\n");
}

void free_list(Node* head) {
    Node* tmp;
    while (head) {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

int main() {
    Node* lista1 = NULL;
    Node* lista2 = NULL;

    srand(time(NULL));

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 0; i < N; i++) {
                int value = rand() % 100;
                #pragma omp task firstprivate(value, i)
                {
                    int tid = omp_get_thread_num();
                    int list_choice = rand() % 2;

                    // 🔒 Região crítica SEM nome — tudo entra aqui
                    #pragma omp critical
                    {
                        if (list_choice == 0) {
                            printf("Thread %d (iteracao %d): INSERINDO %d na Lista 1...\n", tid, i, value);
                            sleep(1); // simula processamento
                            insert(&lista1, value);
                            printf("Thread %d (iteracao %d): FINALIZOU insercao na Lista 1\n", tid, i);
                        } else {
                            printf("Thread %d (iteracao %d): INSERINDO %d na Lista 2...\n", tid, i, value);
                            sleep(1);
                            insert(&lista2, value);
                            printf("Thread %d (iteracao %d): FINALIZOU insercao na Lista 2\n", tid, i);
                        }
                    }
                }
            }
        }
    }

    printf("\nListas finais:\n");
    print_list(lista1, "Lista 1");
    print_list(lista2, "Lista 2");

    free_list(lista1);
    free_list(lista2);

    return 0;
}*/