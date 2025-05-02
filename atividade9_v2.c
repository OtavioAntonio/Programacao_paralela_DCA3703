#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define N 10 // Total de inserções

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
    int num_listas;

    printf("Quantas listas deseja usar? ");
    scanf("%d", &num_listas);

    Node** listas = malloc(num_listas * sizeof(Node*));

    omp_lock_t* locks = malloc(num_listas * sizeof(omp_lock_t));

    for (int i = 0; i < num_listas; i++) {
        listas[i] = NULL;
        omp_init_lock(&locks[i]);
    }

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
                    int list_idx = rand() % num_listas;

                    printf("Thread %d (iteracao %d): INSERINDO %d na Lista %d...\n", tid, i, value, list_idx);
                    
                    omp_set_lock(&locks[list_idx]);
                    insert(&listas[list_idx], value);
                    omp_unset_lock(&locks[list_idx]);

                    printf("Thread %d (iteracao %d): FINALIZOU insercao na Lista %d\n", tid, i, list_idx);
                }
            }
        }
    }

    printf("\nListas finais:\n");
    for (int i = 0; i < num_listas; i++) {
        char name[20];
        sprintf(name, "Lista %d", i);
        print_list(listas[i], name);
        free_list(listas[i]);
        omp_destroy_lock(&locks[i]);
    }

    free(listas);
    free(locks);

    return 0;
}
