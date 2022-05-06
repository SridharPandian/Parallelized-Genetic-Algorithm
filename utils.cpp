#include <float.h>
#include <stdio.h>
#include <string.h>

// Function to return random int between start and end
int get_rand_int(int start, int end) {
    return start + rand() % (end - start); 
}

// Function to return random double between start and end
double get_rand_double(double start, double end) {
    double frac = (double) rand() / (double) RAND_MAX;
    return start + frac*(end - start);  
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

int * get_fittest_individuals(double objective (double *), int num_inds, 
    double ** pop, int pop_size, int pos) {
    
    int * fittest_individuals = (int*) malloc(num_inds * sizeof(int)); 
    int * flags = (int*) calloc(pop_size, sizeof(int)); 

    int min_idx = 0;
    double min_score = DBL_MAX; 
    
    for (int i = 0; i < num_inds; i++) {
        for (int j = 0; j < pop_size; j++) {
            double score = pos * objective(pop[j]); 
            if (score < min_score && flags[j] != -1) {
                min_idx = j; 
                min_score = score; 
            }
        }
        fittest_individuals[i] = min_idx;
        flags[min_idx] = -1; 
        min_score = DBL_MAX; 
    }

    return fittest_individuals; 
}

void swap_fittest_worst (int num_inds, double ** pop_curr, double ** pop_next, 
    int * fittest, int * worst, int num_variables) {
    for (int i = 0; i < num_inds; i++) {
        memcpy(pop_next[worst[i]], pop_curr[fittest[i]], num_variables*sizeof(double));
    }
}

/*
Initializes population of each island for multi-island GA
*/
void initialize_population(double *** pop, double ** bounds, int num_islands, 
    int pop_size, int num_variables) {
    for (int k = 0; k < num_islands; k++) {
        pop[k] = (double**) malloc(pop_size*sizeof(double*));
        #pragma omp parallel for
        for (int i = 0; i < pop_size; i++) {
            pop[k][i] = (double*) malloc(num_variables*sizeof(double));
            for (int j = 0; j < num_variables; j++) {
                pop[k][i][j] = bounds[0][j];
                pop[k][i][j] += get_rand_double(0, 1)*(bounds[1][j] - bounds[0][j]);
            }
        }
    }
}
