#include<iostream>

#include "Individual.h"
#include "Population.h"
#include "GeneticAlgorithm.h"

using namespace std;

int main() {
    GeneticAlgorithm* geneticAlgorithm = new GeneticAlgorithm();

    geneticAlgorithm->population->calcFitness();

    cout << "Generation: " << geneticAlgorithm->generation << " Fittest: " << geneticAlgorithm->population->fittest;

    //While population gets an individual with maximum fitness
    while (geneticAlgorithm->population->fittest < 5) {
        ++geneticAlgorithm->generation;

        //Do selection
        geneticAlgorithm->selection();

        //Do crossover
        geneticAlgorithm->crossOver();

        //Do mutation under a random probability
        if (rand() % 7 < 5) {
            geneticAlgorithm->mutation();
        }

        //Add fittest offspring to population
        geneticAlgorithm->addFittestOffspring();

        //Calculate new fitness value
        geneticAlgorithm->population->calcFitness();

        cout << "Generation: " << geneticAlgorithm->generation << " Fittest: " << geneticAlgorithm->population->fittest;
    }

    cout << "\nSolution found in generation " << geneticAlgorithm->generation;

    // Printing the information about the fittest individual
    int fittestIndividualIdx = geneticAlgorithm->population->fittestIdx;
    Individual fittestIndividual = geneticAlgorithm->population->individuals[fittestIndividualIdx];
    cout << "Fitness: " << fittestIndividual.fitness;

    // Printing out the genes of the fittest individual
    cout << "Genes: ";
    for (int i = 0; i < 5; i++) {
        cout << fittestIndividual.genes[i];
    }
}