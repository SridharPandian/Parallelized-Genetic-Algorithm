#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <random>
#include "utils.cpp"

int perform_selection(int k, int pop_size, double * scores) {
    int rand_idx = get_rand_int(0, pop_size-1);  
    int min_idx = rand_idx; 
    
    for (int i = 0; i < k-1; i++) {
        rand_idx = get_rand_int(0, pop_size-1); 
        min_idx = (scores[rand_idx] < scores[min_idx]) ? rand_idx : min_idx; 
    }

    return min_idx; 
}

void selection(int pop_size, int num_variables, double ** population, double * scores) {
    int k = 3; 
    for (int i = 0; i < pop_size; i++) {
        int next_idx = perform_selection(k, pop_size, scores); 
        // printf("Next idx: %d\n", next_idx); 
        for (int j = 0; j < num_variables; j++)
            population[i][j] = population[next_idx][j]; 
    }
}

void perform_crossover(int num_variables, double * p1, double * p2, double r_cross) {
    if (get_rand_double(0., 1.) < r_cross) {
        std::default_random_engine generator(get_rand_int(0, RAND_MAX));
        std::normal_distribution<double> distribution(0.5,0.15);

        double gamma = distribution(generator); 
        // fprintf(stderr, "Gamma : %f\n", gamma); 
        for (int j = 0; j < num_variables; j++) {
            double temp = p1[j];
            p1[j] = (1 - gamma) * p1[j] + gamma * p2[j]; 
            p2[j] = (1 - gamma) * temp + gamma * p2[j]; 
        }
    }
}

void crossover(int pop_size, int num_variables, double ** population, double r_cross) {
    for (int i = 0; i < pop_size; i += 2) {
        perform_crossover(num_variables, population[i], population[i+1], r_cross); 
    }
}

void mutate(int num_variables, double * p, double ** bounds, double r_mut) {
    if (get_rand_double(0., 1.) < r_mut) { 
        for (int j = 0; j < num_variables; j++) {
            double range = (bounds[1][j] - bounds[0][j]) / 100.0; 
            p[j] += get_rand_double(-1*range, range); 
            
            // limit mutation to bounds
            p[j] = (p[j] > bounds[1][j]) ? bounds[1][j] : p[j]; 
            p[j] = (p[j] < bounds[0][j]) ? bounds[0][j] : p[j];
        }
    }
}

void evaluate_population(double objective (double *), int pop_size, 
    double ** population, double * scores) {
    for (int i = 0; i < pop_size; i++) {
        scores[i] = objective(population[i]); 
    }
}

/*
Genetic algorithm 
*/
void genetic_algorithm(double objective (double *),  double ** bounds, int num_variables, 
    int pop_size, int num_gens, double r_cross, double r_mut) {
    
    if (pop_size % 2 == 1) pop_size += 1; 
    double ** population = (double**) malloc(pop_size*sizeof(double*));
    double * scores = (double*) malloc(pop_size*sizeof(double)); 
    // printf("Before first selection: \n"); 

    srand(time(NULL));
    for (int i = 0; i < pop_size; i++) {
        population[i] = (double*) malloc(num_variables*sizeof(double));
        for (int j = 0; j < num_variables; j++) {
            population[i][j] = bounds[0][j] + get_rand_double(0, 1)*(bounds[1][j] - bounds[0][j]);
            //printf("Lower bound: %f, Upper Bound: %f\n", bounds[0][j], bounds[1][j]);  
            // printf("Pop[i][%d]: %f, ", j, population[i][j]); 
        }
        // printf("\n"); 

        // printf("Score of individual i: %f\n", scores[i]); 
    }

    // printf("After first selection: \n"); 

    for (int i = 0; i < num_gens; i++) {
        evaluate_population(objective, pop_size, population, scores); 
        selection(pop_size, num_variables, population, scores);
        crossover(pop_size, num_variables, population, r_cross); 
        for (int i = 0; i < pop_size; i++)
            mutate(num_variables, population[i], bounds, r_mut);
        
        double min_score = scores[0];
        double * min_individual = population[0]; 
        for (int i = 1; i < pop_size; i++) {
            if (scores[i] < min_score) {
                min_score = scores[i];
                min_individual = population[i];
            }
        }    
        printf("Best score: %f\n", min_score); 
        for (int j = 0; j < num_variables; j++)
            printf("Best individual: Pop[i][%d]: %f, ", j, min_individual[j]); 
        printf("\n");

    }
}