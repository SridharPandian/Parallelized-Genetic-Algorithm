#include "functions.cpp"
#include "multi_island_ga.cpp"

int main() {
    // declare and initialize hyperparameters for GA
    int num_islands, migration_interval, pop_size, num_variables, num_gens; 
    double r_cross, r_mut, r_mig; 
    num_islands = 8; 
    num_variables = 5; 
    pop_size = 1000; 
    num_gens = 100;
    migration_interval = 10; 
    r_cross = 0.8;
    r_mut = 0.15; 
    r_mig = 0.01;

    // get bounds for optimization parameters
    double ** bounds = get_bounds(num_variables, -(10^6), 10^6); 

    // call GA with given hyperparameters
    multi_island_ga_migration(rosenbrock_function_dim5, bounds, 
        num_variables, num_islands, pop_size, num_gens, migration_interval, 
        r_cross, r_mut, r_mig);

    return 0; 
}