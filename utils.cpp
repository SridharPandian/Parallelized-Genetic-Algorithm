#include <float.h>
#include <stdio.h>
#include <string.h>

// Function to return random int between start and end
long int get_rand_int(int start, int end, long int seed) {
    return start + seed % (end - start); 
}

// Function to return random double between start and end
double get_rand_double(double start, double end, long int seed) {
    double frac = (double) seed / (double) RAND_MAX;
    return start + frac*(end - start);  
}

// Function to obtain a probability value using an input seed
double get_prob(long int seed) {
    return (double) seed / (double) RAND_MAX;
}

// Function to change the existing random integer value
long int change_rand_int_value(long int rand_val) {
    return (rand_val * 17931 + 7391) % RAND_MAX;
}

// Function to change the existing random double value
double change_rand_double_value(double rand_val) {
    return (rand_val * 17931 + 7391) / RAND_MAX;
}

// Function to generate a random block containing random values for mutation
long int** get_rand_block(int pop_size, int num_variables) {
    long int** rand_block = (long int**) malloc(pop_size*sizeof(long int*));
    for (int i = 0; i < pop_size; i++) {
        rand_block[i] = (long int*) malloc(num_variables * sizeof(long int));
        
        for (int j = 0; j < num_variables; j++) {
            rand_block[i][j] = rand();
        }
    }
    return rand_block;
}

// Function to generate a random block containing gamma probabilities for crossover
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

// Function to return indexes of numinds fittest individuals of population
int * get_fittest_individuals(double objective (double *), int num_inds, 
    double ** pop, int pop_size, int pos) {

    int * fittest_individuals = (int*) malloc(num_inds * sizeof(int)); 
    int * flags = (int*) calloc(pop_size, sizeof(int)); 

    int min_idx = 0;
    double min_score = DBL_MAX; 
    
    for (int i = 0; i < num_inds; i++) {
        for (int j = 0; j < pop_size; j++) {
            double score = pos * objective(pop[j]); 
            // check if score is below minimum and index not already added
            if (score < min_score && flags[j] != -1) {
                min_idx = j; 
                min_score = score; 
            }
        }
        // add index, set flag and reset min_score
        fittest_individuals[i] = min_idx;
        flags[min_idx] = -1; 
        min_score = DBL_MAX; 
    }

    free(flags);

    return fittest_individuals; 
}

// Function to swap fittest individuals from pop_curr with least fit individuals
// from pop_next 
void swap_fittest_worst (int num_inds, double ** pop_curr, double ** pop_next, 
    int * fittest, int * worst, int num_variables) {
    for (int i = 0; i < num_inds; i++) {
        // do memcpy using fittest and worst indexes
        int f, w = fittest[i], worst[i];
        for (int j = 0; j < num_variables; j++) {
            pop_next[w][j] = pop_curr[f][j]; 
        }  
        //memcpy(pop_next[w], pop_curr[f], num_variables*sizeof(double));
    }
}

// Function to initializes population of each island for multi-island GA
void init_population(double *** pop, double ** bounds, int num_islands, 
    int pop_size, int num_variables) {
    // loop over islands
    for (int k = 0; k < num_islands; k++) {
        pop[k] = (double**) malloc(pop_size*sizeof(double*));
        // loop over population omp
        for (int i = 0; i < pop_size; i++) {
            pop[k][i] = (double*) malloc(num_variables*sizeof(double));
            // randomly initialize between bounds
            for (int j = 0; j < num_variables; j++) {
                pop[k][i][j] = bounds[0][j] + get_prob(rand())*(bounds[1][j] - bounds[0][j]);
            }
        }
    }
}
