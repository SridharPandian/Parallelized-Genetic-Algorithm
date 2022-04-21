#include "genetic_algorithm.cpp"

// test objective funtion
double test_objective(double * w) {
    double offsets[3] = {-8., 17., 22.};
    double objective = 0; 
    for (int i = 0; i < 3; i++)
        objective += pow((w[i] - offsets[i]), 2.); 
    
    return objective;
}

// goldstein-price objective function
double goldstein_price_function(double * w) {
    double x = w[0], y = w[1]; 
    double b1 = pow((x + y + 1), 2); 
    double b2 = 19 - 14*x + 3*pow(x, 2) - 14*y + 6*x*y + 3*pow(y, 2); 
    double t1 = 1 + b1*b2;

    double b3 = pow(2*x - 3*y, 2); 
    double b4 = 18 - 32*x + 12*pow(x, 2) + 48*y - 36*x*y + 27*pow(y, 2); 
    double t2 = 30 + b3*b4;

    return (t1*t2);
}

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
                bounds[i][j] = -2.; 
            else 
                bounds[i][j] = 2.;
        } 
    }

    // call GA with given hyperparameters
    genetic_algorithm(goldstein_price_function, bounds, num_variables, 
        pop_size, num_gens, r_cross, r_mut);     
    
    return 0; 
}