#ifndef OTM_H_
#define OTM_H_
#include <stdio.h>
#include <math.h>

class Optimization{
    public:
        Optimization();
        ~Optimization();

        double function(double x1, double x2);
        double phi(double x1, double x2, double d1, double d2, double t);

        double derivativeX1(double x1, double x2);
        double derivativeX2(double x1, double x2);

        double goldenSection(double x1, double x2, double eps, double p, double d1, double d2);
        double armijo();
        double gradient(double x1, double x2);
        double newton();
        double quasiNewton();
};

#endif