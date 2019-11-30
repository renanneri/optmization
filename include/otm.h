#ifndef OTM_H_
#define OTM_H_
#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;

class Optimization{
    public:
        Optimization();
        ~Optimization();

        double goldenSection(double x1, double x2, double eps, double p, double d1, double d2);
        double armijo(double x1, double x2, double d1, double d2, double gama=0.8, double eta=0.25);
        double gradient(double x1, double x2);
        double newton(double x1, double x2);
        double quasiNewton(double x1, double x2);
};

#endif