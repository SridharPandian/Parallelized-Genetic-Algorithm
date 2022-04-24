#include "genetic_algorithm.cpp"
#include "functions.cpp"

int main() {
    // initialize bounds for GA
    double ** bounds = (double**) malloc(2*sizeof(double*));

    // declare and initialize hyperparameters for GA
    int pop_size, num_variables, num_gens; 
    double r_cross, r_mut; 
    num_variables = 2; 
    pop_size = 250; 
    num_gens = 1000;
    r_cross = 0.8;
    r_mut = 0.1; 

    // initialize bounds for GA
    for (int i = 0; i < 2; i++) {
        bounds[i] = (double*) malloc(num_variables*sizeof(double));
        for (int j = 0; j < num_variables; j++) {
            if (i == 0)
                bounds[i][j] = -100.; 
            else 
                bounds[i][j] = 100.;
        } 
    }

    // call GA with given hyperparameters
    double * opt_res = genetic_algorithm(easom_function, bounds, num_variables, 
        pop_size, num_gens, r_cross, r_mut);

    // print optimized parameters
    for (int j = 0; j < num_variables; j++) {
        printf("Optimized_Result[%d]: %f, ", j, opt_res[j]); 
    }
    printf("\n"); 

    // print minimum score
    printf("Minimum score: %f\n", easom_function(opt_res)); 
    return 0; 
}