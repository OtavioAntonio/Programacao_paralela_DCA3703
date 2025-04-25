#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

typedef struct Node {
    char name[50];
    struct Node* next;
} Node;

int main() {
    const char* stars[] = {
        "Sirius", "Betelgeuse", "Rigel", "Vega", "Aldebaran",
        "Antares", "Procyon", "Altair", "Deneb", "Canopus"
    };
    int num_stars = sizeof(stars) / sizeof(stars[0]);

    // Criar lista encadeada
    Node* head = NULL;
    Node* tail = NULL;
    for (int i = 0; i < num_stars; i++) {
        Node* new_node = (Node*) malloc(sizeof(Node));
        strcpy(new_node->name, stars[i]);
        new_node->next = NULL;

        if (head == NULL) {
            head = new_node;
            tail = new_node;
        } else {
            tail->next = new_node;
            tail = new_node;
        }
    }
    
    #pragma omp parallel
    {
        #pragma omp single
        {
            Node* current = head;
            while (current != NULL) {
                Node* node = current;
                #pragma omp task //firstprivate(node)
                {
                    printf("Arquivo: %s | Thread: %d\n", node->name, omp_get_thread_num());
                }
                current = current->next;
            }
        }
    }
    
    // Liberar memória
    Node* temp;
    while (head != NULL) {
        temp = head;
        head = head->next;
        free(temp);
    }

    return 0;
}