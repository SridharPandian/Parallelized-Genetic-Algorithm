#include <time.h>   
#include <omp.h>
#include "genetic_algorithm.cpp"

/*
Main function to run Multi-Island Genetic algorithm for a given objective function, 
population size, number of generations and other hyperparamters
@param double objective (double *) - objective function to be minimized 
@param double * bounds - bounds for optimization paramters
@param int num_variables - number of optimization paramters
@param int num_islands - number of islands (typically equal to number of threads)
@param int pop_size - size of population
@param int num_gens - number of generations to run for
@param int migration_interval - interval for each migration
@param double r_cross - hyperparameter for crossover rate
@param double r_mut - hyperparameter for mutation rate 
@param double r_mig - hyperparameter for migration rate 
*/
void multi_island_ga(double objective (double *),  double ** bounds, int num_variables,
    int num_islands, int pop_size, int num_gens, int migration_interval, double r_cross, double r_mut,
    double r_mig) {  
    // initialize random seed using time
    srand(time(NULL));

    // always use even population size
    if (pop_size % 2 == 1) pop_size += 1;

    // declare arrays to collect optimized final populations
    double *** pop = (double***) malloc(num_islands*sizeof(double**)); 
    
    // Initializing random blocks
    long int *** sel_rand_block = (long int***) malloc(num_islands*sizeof(long int**)); 
    double ** gamma_block = (double **) malloc(num_islands*sizeof(double*));
    long int*** cross_prob = (long int***) malloc(num_islands*sizeof(long int**));
    long int *** mut_rand_block = (long int***) malloc(num_islands*sizeof(long int**)); 
    long int *** mut_prob = (long int***) malloc(num_islands*sizeof(long int**)); 

    // Assigning all the random values
    for (int i = 0; i < num_islands; i++) {
        // Random values for selection
        sel_rand_block[i] = get_rand_block(pop_size, 3);

        // Random values for crossover
        gamma_block[i] = get_gamma_block(pop_size / 2);
        cross_prob[i] = get_rand_block(pop_size / 2, 1);

        // Random values for mutation
        mut_rand_block[i] = get_rand_block(pop_size, num_variables);
        mut_prob[i] = get_rand_block(pop_size, 1);
    }
    
    // compute total number of migrations
    int num_migrations = num_gens / migration_interval;

    // Initializing the number of threads - incase of serialization
    int thread_num = 1;

    // Initializing the number of threads for the parallel case
    #if defined(_OPENMP)
    omp_set_nested(true);
    printf("Using %d threads!\n", omp_get_max_threads());
    thread_num = omp_get_max_threads() / num_islands;
    printf("Number of threads per island: %d\n", thread_num);
    #endif   

    // set rate of migration and declare arrays to hold fittest and worst
    int num_inds = r_mig * pop_size;  
    int ** fittest = (int**) malloc(num_islands*sizeof(int*));
    int ** worst = (int**) malloc(num_islands*sizeof(int*)); 

    // initialize population of all islands
    init_population(pop, bounds, num_islands, pop_size, num_variables); 

    // Running the genetic island algorithm
    for (int m = 0; m < num_migrations; m++) {
        // call genetic algorithm on all islands with omp parallel
        #pragma omp parallel for num_threads(num_islands)
        for (int i = 0; i < num_islands; i++) {
            pop[i] = genetic_algorithm_multi(objective, pop[i], sel_rand_block[i], gamma_block[i], cross_prob[i], mut_rand_block[i], mut_prob[i], 
                bounds, num_variables, pop_size, migration_interval, r_cross, r_mut, thread_num);

            // get fittest and worst individuals of current population
            fittest[i] = get_fittest_individuals(objective, num_inds, pop[i], pop_size, 1);
            worst[i] = get_fittest_individuals(objective, num_inds, pop[i], pop_size, -1);
        }
        #pragma omp barrier

        // swap fittest of current island with worst of next island - ring topology
        for (int i = 0; i < num_islands; i++) {
            // if i = last island, then next = 0th island
            int next = (i == num_islands-1) ? 0 : i+1;
            swap_fittest_worst(num_inds, pop[i], pop[next], fittest[i], 
                worst[next], num_variables); 
        }
    }

    // Obtaining the fittest score from each island
    #pragma omp parallel for
    for (int i = 0; i < num_islands; i++) {
        // get fittest on island and print objective
        int * fittest_curr = get_fittest_individuals(objective, 1, pop[i], pop_size, 1); 
        int f = fittest_curr[0]; 
        printf("Best score, island %d: %f\n", i, objective(pop[i][f])); 
    }

    // Freeing all the allocated memory
    free(pop);
    free(sel_rand_block);
    free(gamma_block);
    free(cross_prob);
    free(mut_rand_block);
    free(mut_prob);
    free(worst); 
    free(fittest);   
}