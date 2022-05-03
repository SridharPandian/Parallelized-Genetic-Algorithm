#include "functions.cpp"
#include "genetic_algorithm.cpp"

int main() {
    // declare and initialize hyperparameters for GA
    int num_variables, num_gens, thread_num, runs; 
    double r_cross, r_mut; 
    num_variables = 1000; 
    long int pop_size = 100000; 
    num_gens = 1;
    r_cross = 0.8;
    r_mut = 0.15; 
    thread_num = 32;
    runs = 1;

    // get bounds for optimization parameters
    double ** bounds = get_bounds(num_variables, -(10^6), 10^6); 

    // Initializing and starting the timer
    Timer t;
    t.tic();

    // call GA with given hyperparameters
    printf("Starting genetic algorithm!\n");
    double* residuals = genetic_algorithm(rosenbrock_function_dim5, bounds, num_variables, pop_size, num_gens, r_cross, r_mut, thread_num);
    printf("Ending genetic algorithm!\n");

    printf("Time taken by the algorithm: %10f\n", t.toc());

    // print optimized parameters
    // for (int j = 0; j < num_variables; j++) {
    //     printf("\nOpt_Res[%d]: %f, ", j, opt_res[j]); 
    // }
    // printf("\n"); 

    // print minimum score
    printf("Minimum score of run: %f\n", rosenbrock_function_dim5(residuals)); 
    return 0; 

    // free(residuals);
}