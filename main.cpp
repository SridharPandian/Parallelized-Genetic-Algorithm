#include "functions.cpp"
#include "multi_island_ga.cpp"
#include "misc.h"

int main() {
    // declare and initialize hyperparameters for GA
    int num_islands, migration_interval, pop_size, num_variables, num_gens; 
    double r_cross, r_mut, r_mig; 
    num_islands = 8; 
    num_variables = 100; 
    pop_size = 100000; 
    num_gens = 1;
    migration_interval = 1; 
    r_cross = 0.8;
    r_mut = 0.15; 
    r_mig = 0.01;

    // get bounds for optimization parameters
    double ** bounds = get_bounds(num_variables, -(10^6), 10^6); 

    // Initializing a timer
    Timer t;
    t.tic();

    // call GA with given hyperparameters
    multi_island_ga(goldstein_price_function, bounds, 
        num_variables, num_islands, pop_size, num_gens, migration_interval, 
        r_cross, r_mut, r_mig);

    printf("Total time taken by the algorithm: %10f\n", t.toc());

    free(bounds);

    return 0; 
}