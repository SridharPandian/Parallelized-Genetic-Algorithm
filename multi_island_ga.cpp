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
    
    // compute total number of migrations
    int num_migrations = num_gens / migration_interval;

    // set rate of migration and declare arrays to hold fittest and worst
    int num_inds = r_mig * pop_size;  
    int ** fittest = (int**) malloc(num_inds*sizeof(int*));
    int ** worst = (int**) malloc(num_inds*sizeof(int*)); 

    // initialize population of all islands
    init_population(pop, bounds, num_islands, pop_size, num_variables); 
    // printf("initialized %d\n", num_migrations); 

    for (int m = 0; m < num_migrations; m++) {
        // call genetic algorithm on all islands with omp parallel
        #pragma omp parallel for
        for (int i = 0; i < num_islands; i++) {
            // printf("Island GA %d start \n", i); 
            pop[i] = genetic_algorithm_multi(objective, pop[i], bounds, num_variables, pop_size, 
                migration_interval, r_cross, r_mut, rand());

            // get fittest and worst individuals of current population
            fittest[i] = get_fittest_individuals(objective, num_inds, pop[i], pop_size, 1);
            worst[i] = get_fittest_individuals(objective, num_inds, pop[i], pop_size, -1);
            // printf("Island GA %d done \n", i);
        }
        
        // swap fittest of current island with worst of next island - ring topology
        for (int i = 0; i < num_islands; i++) {
            // if i = last island, then next = 0th island
            int next = (i == num_islands-1) ? 0 : i+1;
            swap_fittest_worst(num_inds, pop[i], pop[next], fittest[i], 
                worst[next], num_variables); 
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < num_islands; i++) {
        // get fittest on island and print objective
        int * fittest_curr = get_fittest_individuals(objective, 1, pop[i], pop_size, 1); 
        int f = fittest_curr[0]; 
        printf("Best score, island %d: %f\n", i, objective(pop[i][f])); 
    }    
}