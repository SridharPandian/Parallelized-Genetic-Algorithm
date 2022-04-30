#include "functions.cpp"
#include "genetic_algorithm.cpp"

int main() {
    // declare and initialize hyperparameters for GA
    int num_variables, num_gens, num_threads; 
    double r_cross, r_mut; 
    num_variables = 1000; 
    long int pop_size = 100000; 
    num_gens = 1;
    r_cross = 0.8;
    r_mut = 0.15; 
    num_threads = 32;

    // get bounds for optimization parameters
    double ** bounds = get_bounds(num_variables, -(10^6), 10^6); 

    // Getting the maximum number of threads allowed by the machine
    // #pragma omp parallel
    // num_threads = omp_get_num_threads();

    // Initializing and starting the timer
    Timer t;
    t.tic();

    // call GA with given hyperparameters
    double * opt_res = genetic_algorithm(rosenbrock_function_dim5, bounds, num_variables, pop_size, num_gens, r_cross, r_mut, num_threads);

    // Computing the time
    double time_taken = t.toc();
    printf("Time taken by the algorithm: %10f\n", time_taken);

    // print optimized parameters
    // for (int j = 0; j < num_variables; j++) {
    //     printf("\nOpt_Res[%d]: %f, ", j, opt_res[j]); 
    // }
    // printf("\n"); 

    // print minimum score
    printf("Minimum score: %f\n", rosenbrock_function_dim5(opt_res)); 
    return 0; 
}