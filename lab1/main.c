#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPS 1.0E-20

// f(x) = e^x - 4x^2 - 3x (вариант 5)
double f1(double x) { return exp(x) - 4. * x * x - 3. * x; }
double f1_p(double x) { return exp(x) - 8. * x - 3.; }

// 100 * x * x - 10000 * x - v (вариант 5 => v = 5)
double f2(double x) { return 100 * x * x - 10000 * x - 5; }
double f2_p(double x) { return 200 * (x - 50); }

void newton_method(double (*f)(double), double (*f_p)(double), double x0) {
    double x = x0;
    int iter = 0;
    double error;
    
    printf("\n   > Newton method\n");
    printf("Iter\t x\t\t\t f(x)\t\t\t Error\n");
    
    do {
        if (fabs(f_p(x)) < 1e-15) {
            printf("Err: fp close to zero iter %d\n", iter);
            break;
        }
        
        double x_new = x - f(x) / f_p(x);
        error = fabs(x_new - x);
        
        printf("%d\t\t %.54f\t %.54f\t %.54f\n", iter, x, f(x), error);
        
        x = x_new;
        iter++;
    } while (error > EPS && iter < 100);
    
    printf("Result: x = %.54f, f(x) = %.54f, iter: %d\n", x, f(x), iter);
}

void iteration_method(double (*f)(double), double (*f_p)(double), double x0, double (*phi)(double)) {
    double x = x0;
    int iter = 0;
    double error;
    
    printf("\n   > Iteration Method\n");
    printf("Iter\t x\t\t\t Error\n");
    
    do {
        double x_new = phi(x);
        error = fabs(x_new - x);
        
        printf("%d\t\t %.54f\t %.54f\n", iter, x, error);
        
        x = x_new;
        iter++;
    } while (error > EPS && iter < 100);
    
    printf("Result: x = %.54f, f(x) = %.54f, iter: %d\n", x, f(x), iter);
}

void chord_method(double (*f)(double), double (*f_p)(double), double a, double b) {
    double x = b;
    int iter = 0;
    double error;
    
    printf("\n   > Chord method\n");
    printf("Iter\t x\t\t\t f(x)\t\t\t Error\n");
    
    do {
        double fx = f(x);
        double fa = f(a);
        
        if (fabs(f(x) - f(a)) < 1e-15) {
            printf("Err: fp close to zero iter %d\n", iter);
            break;
        }
        
        double x_new = a - f(a) * (x - a) / (f(x) - f(a));
        error = fabs(x_new - x);
        
        printf("%d\t\t %.54f\t %.54f\t %.54f\n", iter, x, fx, error);
        
        x = x_new;
        iter++;
    } while (error > EPS && iter < 100);
    
    printf("Result: x = %.54f, f(x) = %.54f, iters: %d\n", x, f(x), iter);
}

double phi1(double x) { return x - 0.4 * (exp(x) - 4.0 * x * x - 3.0 * x); }
double phi2(double x) { return (-3.0 + sqrt(9.0 + 16.0 * exp(x))) / 8.0; }
double phi3(double x) { return x - 0.015 * (exp(x) - 4.0 * x * x - 3.0 * x); }

double phi12(double x) { return (100.0 * x * x - 5.0) / 10000.0; }
double phi22(double x) { return sqrt(100.0 * x + 0.05); }

int main() {
    printf("exp1");
    printf("\nRoot 1 = (-0.87031, 0) (interval [-1, 0])");
    newton_method(f1, f1_p, -1);
    printf("\nRoot 2 = (0.3218, 0) (interval [0, 1])");
    newton_method(f1, f1_p, 0);
    printf("\nRoot 3 = (4.58226, 0) (interval [4, 5])");
    newton_method(f1, f1_p, 4);    
    
    printf("\nRoot 1 = (-0.87031, 0) (interval [-1, 0])");
    iteration_method(f1, f1_p, -1, phi1);
    printf("\nRoot 2 = (0.3218, 0) (interval [0, 1])");
    iteration_method(f1, f1_p, 0, phi2);
    printf("\nRoot 3 = (4.58226, 0) (interval [4, 5])");
    iteration_method(f1, f1_p, 4, phi3);
    
    printf("\nRoot 1 = (-0.87031, 0) (interval [-1, 0])");
    chord_method(f1, f1_p, -1, 0);
    printf("\nRoot 2 = (0.3218, 0) (interval [0, 1])");
    chord_method(f1, f1_p, 0, 1);
    printf("\nRoot 3 = (4.58226, 0) (interval [4, 5])");
    chord_method(f1, f1_p, 4, 5);

    printf("\n\nexp3");
    printf("\nRoot 1 = (-0.0005, 0) (interval [-1, 0])");
    newton_method(f2, f2_p, -1);
    iteration_method(f2, f2_p, -1, phi12);
    chord_method(f2, f2_p, -1, 0);

    printf("\nRoot 2 = (100.0005, 0) (interval [100, 101])");
    newton_method(f2, f2_p, 100);
    iteration_method(f2, f2_p, 100, phi22);
    chord_method(f2, f2_p, 100, 101);
}