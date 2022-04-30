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