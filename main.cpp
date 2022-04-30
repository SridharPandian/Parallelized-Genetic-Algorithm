#include "functions.cpp"
#include "genetic_algorithm.cpp"

int main() {
    // declare and initialize hyperparameters for GA
    int pop_size, num_variables, num_gens; 
    double r_cross, r_mut; 
    num_variables = 5; 
    pop_size = 10000; 
    num_gens = 1000;
    r_cross = 0.8;
    r_mut = 0.15; 

    // get bounds for optimization parameters
    double ** bounds = get_bounds(num_variables, -(10^6), 10^6); 

    // call GA with given hyperparameters
    double * opt_res = genetic_algorithm(rosenbrock_function_dim5, 
        bounds, num_variables, pop_size, num_gens, r_cross, r_mut);

    // print optimized parameters
    for (int j = 0; j < num_variables; j++) {
        printf("Opt_Res[%d]: %f, ", j, opt_res[j]); 
    }
    printf("\n"); 

    // print minimum score
    printf("Minimum score: %f\n", rosenbrock_function_dim5(opt_res)); 
    return 0; 
}