#include <stdlib.h>

#define GENE_LENGTH 5

class Individual {
    public:
        // Gene information
        int geneLength = GENE_LENGTH;
        int* genes = (int*) malloc(geneLength * sizeof(int));

        // Fitness value
        int fitness = 0;

        // Constructor function
        Individual() {
            for (int idx = 0; idx < geneLength; idx++)
                genes[idx] = rand() % 2;

            fitness = 0;
        }

        // Function to calculate the fitness
        void calcFitness() {
            fitness = 0;

            // Fitness value is the number of genes which are 1
            for (int idx = 0; idx < geneLength; idx++) {
                if (genes[idx] == 1) {
                    ++fitness;
                }
            }
        }
};