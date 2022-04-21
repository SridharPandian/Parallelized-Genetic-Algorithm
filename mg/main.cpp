#include "genetic_algorithm.cpp"
#include "functions.cpp"

int main() {
    // initialize bounds for GA
    double ** bounds = (double**) malloc(2*sizeof(double*));

    // declare and initialize hyperparameters for GA
    int pop_size, num_variables, num_gens; 
    double r_cross, r_mut; 
    num_variables = 2; 
    pop_size = 100; 
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
    genetic_algorithm(easom_function, bounds, num_variables, 
        pop_size, num_gens, r_cross, r_mut);     
    
    return 0; 
}