#include <stdlib.h>
#include "Population.h"
#include "Individual.h"

#define POP_SIZE 10

class GeneticAlgorithm {
    public:
        int generation = 0;
        Population* population;
        Individual fittestIndividual;
        Individual secondFittestIndividual;

        GeneticAlgorithm () {
            population = new Population(POP_SIZE);
        }

        void selection() {
            // Selecting the most fittest individual
            int fittestIndividualIdx = population->getFittestIndividualIdx();
            fittestIndividual = population->getIndividual(fittestIndividualIdx);

            // Selecting the second most fittest individual
            int fittestIndividualIdx = population->getSecondFittestIndividualIdx();
            secondFittestIndividual = population->getIndividual(fittestIndividualIdx);
        }

        void crossOver() {
            // Selecting a random crossover point
            int crossOverPoint = rand() % population->getIndividual(0).geneLength;

            // Swapping gene values among parents
            for (int idx = 0; idx < crossOverPoint; idx++) {
                int swapValue = fittestIndividual.genes[idx];
                fittestIndividual.genes[idx] = secondFittestIndividual.genes[idx];
                secondFittestIndividual.genes[idx] = swapValue;
            }
        }

        void mutation() {
            int MutationPoint = rand() % population->getIndividual(0).geneLength;

            // Flipping values at the mutation point in the fittest individual
            if (fittestIndividual.genes[MutationPoint] == 0) {
                fittestIndividual.genes[MutationPoint] = 1;
            } else {
                fittestIndividual.genes[MutationPoint] = 0;
            }

            int MutationPoint = rand() % population->getIndividual(0).geneLength;

            // Flipping values at the mutation point in the second fittest individual
            if (secondFittestIndividual.genes[MutationPoint]) {
                secondFittestIndividual.genes[MutationPoint] = 1;
            } else {
                secondFittestIndividual.genes[MutationPoint] = 0;
            }
        }

        Individual getFittestOffspring() {
            if (fittestIndividual.fitness > secondFittestIndividual.fitness)
                return fittestIndividual;
            return secondFittestIndividual; 
        }

        void addFittestOffspring() {
            fittestIndividual.calcFitness();
            secondFittestIndividual.calcFitness();

            int leastFittestIndividualIdx = population->getLeastFittestIndividualIdx();

            population->individuals[leastFittestIndividualIdx] = getFittestOffspring();
        }
};