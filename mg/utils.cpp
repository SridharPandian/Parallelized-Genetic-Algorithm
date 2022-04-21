int get_rand_int(int start, int end) {
    return start + rand() % (end - start); 
}

double get_rand_double(double start, double end) {
    double frac = (double) rand() / (double) RAND_MAX;
    return start + frac*(end - start);  
}