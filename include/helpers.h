#ifndef HELPERS_H_
#define HELPERS_H_
#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;

double function(double x1, double x2);
double phi(double x1, double x2, double d1, double d2, double t);

double derivativeX1(double x1, double x2);
double derivativeX2(double x1, double x2);

double hessianX1X1(double x1, double x2);
double hessianX1X2(double x1, double x2);
double hessianX2X1(double x1, double x2);
double hessianX2X2(double x1, double x2);

vector<double> dnewton(double x1, double x2,double gradX1, double gradX2);
vector<vector<double> > BFGS(vector<vector<double> > H, double x1,double x2,double x1_1, double x2_1);

vector<vector<double>> inverse(vector<vector<double> > H);

vector<double> p(double x1,double x2,double x1_1, double x2_1);
vector<double> q(double x1,double x2,double x1_1, double x2_1);

#endif