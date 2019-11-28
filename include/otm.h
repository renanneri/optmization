#ifndef OTM_H_
#define OTM_H_
#include <stdio.h>
#include <math.h>
#include <vector>

class Optimization{
    public:
        Optimization();
        ~Optimization();

        double function(double x1, double x2);
        double phi(double x1, double x2, double d1, double d2, double t);

        double derivativeX1(double x1, double x2);
        double derivativeX2(double x1, double x2);

        double hessianX1X1(double x1, double x2);
        double hessianX1X2(double x1, double x2);
        double hessianX2X1(double x1, double x2);
        double hessianX2X2(double x1, double x2);

        double goldenSection(double x1, double x2, double eps, double p, double d1, double d2);
        double armijo();
        double gradient(double x1, double x2);
        double newton(double x1, double x2);
        vector<double> dnewton(double x1, double x2,double gradX1, double gradX2);
        double quaseNewton(double x1, double x2);
        vector<vector<double>> BFGS(vector<vector<double>> H, double x1,double x2,double x1_1, double x2_1);

        vector<double> p(double x1,double x2,double x1_1, double x2_1);
        vector<double> q(double x1,double x2,double x1_1, double x2_1);
};

#endif