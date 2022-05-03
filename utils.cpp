#include <random>
#include <math.h>

// Function to return random int between start and end
int get_rand_int(int start, int end, int seed) {
    return start + seed % (end - start); 
}

// Function to return random double between start and end
double get_rand_double(double start, double end, int seed) {
    double frac = (double) seed / (double) RAND_MAX;
    return start + frac*(end - start);  
}

double get_prob(int seed) {
    return (double) seed / (double) RAND_MAX;
}

int** get_rand_block(int pop_size, int num_variables) {
    int** rand_block = (int**) malloc(pop_size*sizeof(int*));
    for (int i = 0; i < pop_size; i++) {
        rand_block[i] = (int*) malloc(num_variables * sizeof(int));
        
        for (int j = 0; j < num_variables; j++) {
            rand_block[i][j] = rand();
        }
    }
    return rand_block;
}

double* get_gamma_block(int block_length) {
    double* gamma_block = (double*) malloc(block_length * sizeof(double));
    for (int i = 0; i < block_length; i++) {
        std::default_random_engine generator(get_rand_int(0, RAND_MAX, rand()));
        std::normal_distribution<double> distribution(0.5,0.15);
        gamma_block[i] = distribution(generator);
    }

    return gamma_block;
}

// Function to generate array of bounds for genetic algorithm 
double ** get_bounds (int num_variables, double lower, double upper) {
    double ** bounds = (double**) malloc(2*sizeof(double*));

    for (int i = 0; i < 2; i++) {
        bounds[i] = (double*) malloc(num_variables*sizeof(double));
        for (int j = 0; j < num_variables; j++) {
            if (i == 0)
                bounds[i][j] = lower; 
            else 
                bounds[i][j] = upper;
        } 
    }

    return bounds; 
}