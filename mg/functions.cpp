#include <math.h>

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

// easom objective function 
double easom_function(double * w) {
    double x = w[0], y = w[1]; 
    double t1 = -1 * cos(x) * cos(y); 

    double b1 = pow(x - M_PI, 2) + pow(y - M_PI, 2); 
    double t2 = exp(-1 * b1); 

    return (t1*t2); 
}