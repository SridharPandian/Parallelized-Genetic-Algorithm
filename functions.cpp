#include <math.h>

// Test objective funtion
// Bounds: -25 <= x, y <= 25
// Min parameters (-8, 17, 22), score: 0
double test_objective(double * w) {
    double offsets[3] = {-8., 17., 22.};
    double objective = 0; 
    for (int i = 0; i < 3; i++)
        objective += pow((w[i] - offsets[i]), 2.); 
    
    return objective;
}

// 5 dimensional Rosenbrock function 
// Bounds: -inf < x_i < inf
// Min parameters: (1, 1, 1, 1, 1), Min score: 0
double rosenbrock_function_dim5(double * w) {
    int n = 5; 
    double sum = 0; 
    for (int i = 0; i < n-1; i++) {
        sum += 100 * pow((w[i+1] - pow(w[i], 2)), 2); 
        sum += pow(1 - pow(w[i], 2), 2); 
    }
    return sum; 
}

// Goldstein-price objective function
// Bounds: -2 <= x, y <= 2
// Min parameters (0, -1), score: 3
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

// Easom objective function 
// Bounds: -100 <= x, y <= 100
// Min parameters (PI, PI), score: -1
double easom_function(double * w) {
    double x = w[0], y = w[1]; 
    double t1 = -1 * cos(x) * cos(y); 

    double b1 = pow(x - M_PI, 2) + pow(y - M_PI, 2); 
    double t2 = exp(-1 * b1); 

    return (t1*t2); 
}

// Holder table objective function
// Bounds: -10 <= x, y <= 10
// Min parameters (+/-8.05502, +/-9.66459)/, Min score: -19.2085
double holder_table_function(double * w) {
    double x = w[0], y = w[1]; 
    double t1 = sin(x) * cos(y); 
    double t2 = abs(1. - sqrt(pow(x, 2) + pow(y, 2))/M_PI);

    return (-1*abs(t1*exp(t2)));  
}

// Eggholder objective function
// Bounds: -512 <= x, y <= 512
// Min parameters (512, 404.2319), Min Score: -959.6407
double eggholder_function(double * w) {
    double x = w[0], y = w[1]; 
    double b1 = sqrt(abs(x/2 + y + 47));
    double t1 = -1 * (y+47) * sin (b1);

    double b2 = sqrt(abs(x - (y + 47)));
    double t2 = -1 * (x) * sin(b2); 

    return (t1+t2); 
}
