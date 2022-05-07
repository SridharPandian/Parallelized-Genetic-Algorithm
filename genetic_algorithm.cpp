#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <random>
#include <limits>
#include "utils.cpp"
#include "misc.h"

/*
Performs a randomization step. All the integer elements of the random block are changed
@param int ** rand_block - random block
@param int block_rows - number of rows in the random block
@param int block_cols - number of columns in the random block
*/
void change_rand_block(long int ** rand_block, int block_rows, int block_cols, int thread_num) {
    #pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < block_rows; i++)
        for (int j = 0; j < block_cols; j++)
            rand_block[i][j] = change_rand_int_value(rand_block[i][j]);
}

/*
Performs a randomization step. All the double elements of the gamma block are changed
@param int ** rand_block - random block
@param int block_rows - number of rows in the random block
@param int block_cols - number of columns in the random block
*/
void change_gamma_block(double * gamma_block, int num_values, int thread_num) {
    #pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < num_values; i++)
        gamma_block[i] = change_rand_double_value(gamma_block[i]);
}

/*
@param int num_variables - number of optimization paramters
@param int pop_size - size of population
@param long int ** sel_rand_block - 2D array containing random values for the selection process
@param double * gamma_block - array containing random crossover gamma probabilities 
@param long int ** cross_prob - 2D array of crossover probabilities
@param long int ** mut_rand_block - 2D array individual gene mutation probabilities
@param long int ** mut_prob - array of individual mutation probabilities
@param int thread_num - number of threads used for parallelization
*/
void rerandomize(int pop_size, int num_variables, long int ** sel_rand_block, double * gamma_block, long int ** cross_prob, long int ** mut_rand_block, long int ** mut_prob, int thread_num) {
    // Change the random values in the random blocks
    change_rand_block(sel_rand_block, pop_size, 3, thread_num);
    change_gamma_block(gamma_block, pop_size / 2, thread_num);
    change_rand_block(cross_prob, pop_size / 2, 1, thread_num);
    change_rand_block(mut_rand_block, pop_size, num_variables, thread_num);
    change_rand_block(mut_prob, pop_size, 1, thread_num);
}

/*
Performs one selection step
@param int k - hyperparamter for selection 
@param int pop_size - size of population
@param double * scores - array of population scores
@param rand_row - array of random values
@return int min_idx - index of individual with 
    minimum score after k selections
*/
int perform_selection(int k, int pop_size, double * scores, long int* rand_row) {
    long int rand_idx = get_rand_int(0, pop_size-1, rand_row[0]);  
    long int min_idx = rand_idx; 

    for (int i = 0; i < k-1; i++) {
        rand_idx = get_rand_int(0, pop_size-1, rand_row[i + 1]); 
        min_idx = (scores[rand_idx] < scores[min_idx]) ? rand_idx : min_idx; 
    }
    return min_idx; 
}

/*
Performs selection for entire population - calls perform_selection()
@param int pop_size - size of population
@param int num_variables - number of optimization paramters
@param double ** pop_curr - array of individuals for population 
                            in current generation
@param double ** pop_curr - array of individuals for population
                            in next generation
@param double * scores - array of population scores
@param int ** rand_block - 2D array of random values
*/
void selection(int pop_size, int num_variables, double ** pop_curr, double ** pop_next, double * scores, long int** rand_block) {
    int k = 3; 
    for (int i = 0; i < pop_size; i++) {
        long int next_idx = perform_selection(k, pop_size, scores, rand_block[i]); 
        for (int j = 0; j < num_variables; j++)
            pop_next[i][j] = pop_curr[next_idx][j];
    }
}

/*
Performs one crossover step to transform two parents into two children
in-place (uses same array as parent to store child)
@param int num_variables - number of optimization paramters
@param double * p1 - parent 1 from current generation
@param double * p2 - parent 2 from current generation
@param double r_cross - hyperparameter for crossover rate
@param double gamma - array containing ramdom gamma values for crossover
@param double cross_prob - random crossover probabilty
*/
void perform_crossover(int num_variables, double * p1, double * p2, double r_cross, double gamma, double cross_prob) {
    if (cross_prob < r_cross) {
        for (int j = 0; j < num_variables; j++) {
            double temp = p1[j];
            p1[j] = (1 - gamma) * p1[j] + gamma * p2[j]; 
            p2[j] = (1 - gamma) * temp + gamma * p2[j]; 
        }
    }
}

/*
Performs crossover for entire population - calls perform_crossover()
@param int pop_size - size of population
@param double ** pop_curr - array of individuals for population 
                            in current generation
@param int num_variables - number of optimization paramters
@param double r_cross - hyperparameter for crossover rate
@param int thread_num - number of threads used for parallelization
@param double * gamma_block - array containing random gamma values
@param long int ** cross_prob - 2D array containing crossover probabilities
*/
void crossover(int pop_size, int num_variables, double ** pop_curr, double r_cross, int thread_num, double* gamma_block, long int** cross_prob) {
    #pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < pop_size; i += 2)
        perform_crossover(num_variables, pop_curr[i], pop_curr[i+1], r_cross, gamma_block[i / 2], get_prob(cross_prob[i / 2][0])); 
}

/*
Performs mutation on a single individual in-place
@param int num_variables - number of optimization paramters
@param double * p - individual from current generation to be mutated
@param double * bounds - bounds for optimization paramters
@param double r_mut - hyperparameter for mutation rate
@param int * rand_row - array of random integers
@param double mut_prob - random probability of mutating an individual
*/
void perform_mutation(int num_variables, double * p, double ** bounds, double r_mut, long int* rand_row, double mut_prob) {
    if (mut_prob < r_mut) { 
        for (int j = 0; j < num_variables; j++) {
            double range = (bounds[1][j] - bounds[0][j]) / 100.0; 
            p[j] += get_rand_double(-1 * range, range, rand_row[j]); 
            
            // limit mutation to bounds
            p[j] = (p[j] > bounds[1][j]) ? bounds[1][j] : p[j]; 
            p[j] = (p[j] < bounds[0][j]) ? bounds[0][j] : p[j];
        }
    }
}

/*
Performs mutation for entire population - calls perform_mutation()
@param int pop_size - size of population
@param double ** pop_curr - array of individuals for population 
                            in current generation
@param int num_variables - number of optimization paramters
@param bounds - bounds for optimization paramters
@param double r_mut - hyperparameter for mutation rate
@param rand_block - 2D array containing random integers
@param int thread_num - number of threads used for parallelization
@param long int ** mut_prob - 2D array containing the mutation probabilities
*/
void mutate(int pop_size, int num_variables, double ** pop_curr, double ** bounds, double r_mut, long int** rand_block, int thread_num, long int** mut_prob) {
    #pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < pop_size; i++)
        perform_mutation(num_variables, pop_curr[i], bounds, r_mut, rand_block[i], get_prob(mut_prob[i][0]));
}

/*
Evaluates scores of current population and return individual with 
minimum score
@param double objective (double *) - objective function to be minimized 
@param int pop_size - size of population
@param double ** pop_curr - array of individuals for population 
                            in current generation
@param double * scores - array of population scores
@param int thread_num - number of threads used for parallelization
@return double * min_individual - individual with minimum score 
*/
double * evaluate_scores(double objective (double *), int pop_size, double ** pop_curr, double * scores, int thread_num) {
    // Initialize score and individual
    double min_score = std::numeric_limits<float>::infinity();
    double * min_individual = pop_curr[0]; 
    
    // Loop over population to get all scores
    #pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < pop_size; i++)
        scores[i] = objective(pop_curr[i]); 

    // Loop over population to get the minimum valued individual
    for (int i = 0; i < pop_size; i++) {   
        if (scores[i] < min_score) {
            min_score = scores[i];
            min_individual = pop_curr[i];
        }
    }

    // Return fittest individual
    return min_individual; 
}

/*
Main function to run Genetic algorithm for a given objective function, 
population size, number of generations and other hyperparamters
@param double objective (double *) - objective function to be minimized 
@param double ** pop_curr - 2D array containing the population values
@param long int ** sel_rand_block - 2D array containing random values for the selection process
@param double * gamma_block - array containing random crossover gamma probabilities 
@param long int ** cross_prob - 2D array of crossover probabilities
@param long int ** mut_rand_block - 2D array individual gene mutation probabilities
@param long int ** mut_prob - array of individual mutation probabilities
@param double * bounds - bounds for optimization paramters
@param int num_variables - number of optimization paramters
@param int pop_size - size of population
@param int num_gens - number of generations to run for
@param double r_cross - hyperparameter for crossover rate
@param double r_mut - hyperparameter for mutation rate 
@param thread_num - number of threads used to parallelize this instance of the algorithm
@return double * pop_curr - final individual with minimum score 
*/
double ** genetic_algorithm_multi(double objective (double *), double ** pop_curr, long int ** sel_rand_block, 
    double * gamma_block, long int ** cross_prob, long int ** mut_rand_block, long int ** mut_prob,
    double ** bounds, int num_variables, int pop_size, int num_gens, double r_cross, double r_mut, int thread_num) { 
    // OMP parameters
    // #if defined(_OPENMP)
    // omp_set_num_threads(thread_num)
    // #endif

    Timer t;

    // Allocate arrays for populations and scores
    double ** pop_next = (double**) malloc(pop_size*sizeof(double*)); 
    double * scores = (double*) malloc(pop_size*sizeof(double));     

    // Initialize population with random values between bounds
    // printf("Initializing the parameters in the memory!\n");
    for (int i = 0; i < pop_size; i++) {
        pop_next[i] = (double*) malloc(num_variables*sizeof(double)); 
    }
    // Evaluate initial scores and get fittest individual
    double * min_individual = evaluate_scores(objective, pop_size, pop_curr, scores, thread_num); 

    // Main GA for loop
    for (int i = 0; i < num_gens; i++) {
        // Changing the random values in all the random blocks
        // t.tic();
        rerandomize(pop_size, num_variables, sel_rand_block, gamma_block, cross_prob, mut_rand_block, mut_prob, thread_num);
        // printf("Time taken for randomization: %10f\n", t.toc());

        
        // Select individuals for next generation
        // t.tic();
        selection(pop_size, num_variables, pop_curr, pop_next, scores, sel_rand_block);
        // printf("Time taken for selection: %10f\n", t.toc());
        
        // Swap next generation with current generation
        std::swap(pop_curr, pop_next); 

        // Perform crossover 
        // t.tic();
        crossover(pop_size, num_variables, pop_curr, r_cross, thread_num, gamma_block, cross_prob); 
        // printf("Time taken for crossover: %10f\n", t.toc());

        // Perform mutation
        // t.tic();
        mutate(pop_size, num_variables, pop_curr, bounds, r_mut, mut_rand_block, thread_num, mut_prob);
        // printf("Time taken for mutation: %10f\n", t.toc());
        
        // Evaluate fitness for population and get fittest individual
        t.tic();
        min_individual = evaluate_scores(objective, pop_size, pop_curr, scores, thread_num); 
        printf("Time taken for evaluation: %10f\n", t.toc());
    }

    // Freeing all the allocated memory
    free(scores);
    free(pop_next);
    free(min_individual);

    // Return final fittest individual
    return pop_curr; 
}