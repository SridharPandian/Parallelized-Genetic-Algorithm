#include <vector>
#include <math.h>

#include "Individual.h"

class Population {
    public:
        // Assigning the population size value
        int popSize;
        
        // Creating an array of Inidividuals
        std::vector<Individual> individuals;

        // Indicator for the fittest individuals
        int fittest = 0;
        int fittestIdx;

        // Constructor function
        Population (int populationSize) {
            popSize = populationSize;

            // Initializing individuals and adding them to the vector array
            for (int idx = 0; idx < popSize; idx++) 
                individuals.push_back(Individual());
        }

        // Obtaining the index of the fittest Individual
        int getFittestIndividualIdx() {
            double maxFit = - INFINITY;
            int maxFitIdx = 0;

            for (int idx = 0; idx < popSize; idx++) {
                if (maxFit <= individuals[idx].fitness) {
                    maxFit = individuals[idx].fitness;
                    maxFitIdx = idx;
                }
            }

            fittest = individuals[maxFitIdx].fitness;
            return maxFitIdx;
        }

        // Obtaining the index of the second fittest individual
        int getSecondFittestIndividualIdx() {
            int maxFit1 = 0;
            int maxFit2 = 0;

            for (int idx = 0; idx < popSize; idx++) {
                if (individuals[idx].fitness > individuals[maxFit1].fitness) {
                    maxFit2 = maxFit1;
                    maxFit1 = idx;
                } else if (individuals[idx].fitness > individuals[maxFit2].fitness) {
                    maxFit2 = idx;
                }
            }

            return maxFit2;
        }

        // Obtaining the index of the least fit individual
        int getLeastFittestIndividualIdx() {
            int minFit = INFINITY;
            int minFitIdx = 0;

            for (int idx = 0; idx < popSize; idx++) {
                if (minFit >= individuals[idx].fitness) {
                    minFit = individuals[idx].fitness;
                    minFitIdx = idx;
                }
            }

            return minFitIdx;
        }

        // Calculate the fitness for all individuals
        void calcFitness() {
            for (int idx = 0; idx < popSize; idx++) {
                individuals[idx].calcFitness();
            }
            fittestIdx = getFittestIndividualIdx();
        }

        Individual getIndividual(int idx) {
            return individuals[idx];
        }
};