#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <random>
#include <limits>
#include "utils.cpp"
#include "misc.h"
#include "omp.h"

/*
Performs one selection step
@param int k - hyperparamter for selection 
@param int pop_size - size of population
@param double * scores - array of population scores
@return int min_idx - index of individual with 
    minimum score after k selections
*/
int perform_selection(int k, int pop_size, double * scores, int* rand_row) {
    int rand_idx = get_rand_int(0, pop_size-1, rand_row[0]);  
    int min_idx = rand_idx; 

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
*/
void selection(int pop_size, int num_variables, double ** pop_curr, double ** pop_next, double * scores) {
    int k = 3; 
    int** rand_block = get_rand_block(pop_size, 3);

    for (int i = 0; i < pop_size; i++) {
        int next_idx = perform_selection(k, pop_size, scores, rand_block[i]); 
        for (int j = 0; j < num_variables; j++)
            pop_next[i][j] = pop_curr[next_idx][j]; 
    }

    free(rand_block);
}

/*
Performs one crossover step to transform two parents into two children
in-place (uses same array as parent to store child)
@param int num_variables - number of optimization paramters
@param double * p1 - parent 1 from current generation
@param double * p2 - parent 2 from current generation
@param double r_cross - hyperparameter for crossover rate
*/
void perform_crossover(int num_variables, double * p1, double * p2, double r_cross, double gamma, int cross_prob) {
    if (get_rand_double(0., 1., cross_prob) < r_cross) {
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
*/
void crossover(int pop_size, int num_variables, double ** pop_curr, double r_cross) {
    double* gamma_block = get_gamma_block(pop_size / 2);
    int** cross_prob = get_rand_block(pop_size / 2, 1);

    #pragma omp parallel for 
    for (int i = 0; i < pop_size; i += 2) {
        perform_crossover(num_variables, pop_curr[i], pop_curr[i+1], r_cross, gamma_block[i / 2], cross_prob[i / 2][0]); 
    }

    free(gamma_block);
    free(cross_prob);
}

/*
Performs mutation on a single individual in-place
@param int num_variables - number of optimization paramters
@param double * p - individual from current generation to be mutated
@param double * bounds - bounds for optimization paramters
@param double r_mut - hyperparameter for mutation rate
*/
void perform_mutation(int num_variables, double * p, double ** bounds, double r_mut, int* rand_row, int mut_prob) {
    if (get_prob(mut_prob) < r_mut) { 
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
*/
void mutate(int pop_size, int num_variables, double ** pop_curr, double ** bounds, double r_mut) {
    Timer t;
    t.tic();
    int** rand_block = get_rand_block(pop_size, num_variables);
    int** mut_prob = get_rand_block(pop_size, 1);
    printf("Time taken to compute random variables: %10f\n", t.toc());

    #pragma omp parallel for 
    for (int i = 0; i < pop_size; i++)
        perform_mutation(num_variables, pop_curr[i], bounds, r_mut, rand_block[i], mut_prob[i][0]);

    free(mut_prob);
    free(rand_block);
}

/*
Evaluates scores of current population and return individual with 
minimum score
@param double objective (double *) - objective function to be minimized 
@param int pop_size - size of population
@param double ** pop_curr - array of individuals for population 
                            in current generation
@param double * scores - array of population scores
@return double * min_individual - individual with minimum score 
*/
double * evaluate_scores(double objective (double *), int pop_size, double ** pop_curr, double * scores) {
    // initialize score and individual
    double min_score = std::numeric_limits<float>::infinity();
    double * min_individual = pop_curr[0]; 
    
    // loop over population
    #pragma omp parallel for
    for (int i = 0; i < pop_size; i++)
        scores[i] = objective(pop_curr[i]); 

    // update minimum individual
    for (int i = 0; i < pop_size; i++) {   
        if (scores[i] < min_score) {
            min_score = scores[i];
            min_individual = pop_curr[i];
        }
    }

    // return fittest individual
    return min_individual; 
}

/*
Main function to run Genetic algorithm for a given objective function, 
population size, number of generations and other hyperparamters
@param double objective (double *) - objective function to be minimized 
@param double * bounds - bounds for optimization paramters
@param int num_variables - number of optimization paramters
@param int pop_size - size of population
@param int num_gens - number of generations to run for
@param double r_cross - hyperparameter for crossover rate
@param double r_mut - hyperparameter for mutation rate 
@return double * min_individual - final individual with minimum score 
*/
double * genetic_algorithm(double objective (double *),  double ** bounds, int num_variables, 
    int pop_size, int num_gens, double r_cross, double r_mut, int thread_num) { 
    // Setting up OpenMP params
    // #if defined(_OPENMP)
    // printf("Using %d threads!\n", thread_num);
    // omp_set_num_threads(thread_num);
    // #endif

    // Initializing timer
    Timer t;

    // always use even population size
    if (pop_size % 2 == 1) pop_size += 1; 
    
    // allocate arrays for populations and scores
    double ** pop_curr = (double**) malloc(pop_size*sizeof(double*));
    double ** pop_next = (double**) malloc(pop_size*sizeof(double*)); 
    double * scores = (double*) malloc(pop_size*sizeof(double)); 

    // initialize random seed
    srand(time(NULL));

    t.tic();
    // initialize population with random values between bounds
    printf("Initializing the parameters in the memory!\n");
    for (int i = 0; i < pop_size; i++) {
        pop_curr[i] = (double*) malloc(num_variables*sizeof(double));
        pop_next[i] = (double*) malloc(num_variables*sizeof(double)); 

        for (int j = 0; j < num_variables; j++)
            pop_curr[i][j] = bounds[0][j] + get_rand_double(0, 1, rand()) * (bounds[1][j] - bounds[0][j]); 
    }
    printf("Time taken for allocating memory: %10f\n", t.toc());

    // evaluate initial scores and get fittest individual
    t.tic();
    printf("Computing the scores for all the individuals in the population before starting the algorithm!\n");
    double * min_individual = evaluate_scores(objective, pop_size, pop_curr, scores); 
    printf("Time taken for first evaluation: %10f\n\n", t.toc());

    // Main GA for loop
    for (int i = 0; i < num_gens; i++) {
        printf("Working on generation: %d\n", i);

        // select individuals for next generation
        t.tic();
        printf("Performing selection for generation: %d\n", i);
        selection(pop_size, num_variables, pop_curr, pop_next, scores);
        printf("Time taken for selection: %10f\n\n", t.toc());
        
        // swap next generation with current generation
        std::swap(pop_curr, pop_next); 

        // perform crossover and mutation
        t.tic();
        printf("Performing crossover for generation: %d\n", i);
        crossover(pop_size, num_variables, pop_curr, r_cross); 
        printf("Time taken for crossover: %10f\n\n", t.toc());

        t.tic();
        printf("Performing mutation for generation: %d\n", i);
        mutate(pop_size, num_variables, pop_curr, bounds, r_mut);
        printf("Time taken for mutation: %10f\n\n", t.toc());
        
        // evaluate fitness for population and get fittest individual
        t.tic();
        printf("Evaluating the scores of each individual after generation: %d!\n", i);
        min_individual = evaluate_scores(objective, pop_size, pop_curr, scores); 
        printf("Time taken for evaluation for this generation: %10f\n\n", t.toc());

    }

    // return final fittest individual
    return min_individual; 
}