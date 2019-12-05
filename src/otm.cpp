#include <stdio.h>
#include <iostream>
#include <math.h> 
#include <string>
#include <vector>
#include "otm.h"
#include "helpers.h"

using namespace std;

Optimization::Optimization(){

}

Optimization::~Optimization(){

}

double Optimization::goldenSection(double x1, double x2, double eps, double p, double d1, double d2)
{
  double theta1 = (3 - sqrt(5))/2.0;
  double theta2 = 1 - theta1;

  // Obtenção do intervalo [a, b]
  double a = 0.0;
  double s = p;
  double b = 2*p;

  // cout << "func " << derivativeX1(x1,x2)*d1 + derivativeX2(x1, x2)*d2 << endl;


  while (function(x1+b*d1, x2+b*d2) < function(x1+s*d1, x2+s*d2)){
    a = s;
    s = b;
    b = 2*b;
    // cout << "f: " << function(x1, x2) << endl;
    // cout << "phi1: " << phi(x1, x2, d1, d2, b) << endl;
    // cout << "phi2: " << phi(x1, x2, d1, d2, s) << endl;
    // cout << "a:" << a << " s:" << s << " b:" << b << endl;
  }

  // Obtenção de t*
  double u = a + theta1*(b - a);
  double v = a + theta2*(b - a);
  int cont = 0;
  while (abs(b - a) > eps){
    // cout << "b - a: " << b-a << endl;
    // if (cont > 1000) break;
    if (function(x1+u*d1, x2+u*d2) < function(x1+v*d1, x2+v*d2)){
      b = v;
      v = u;
      u = a + theta1*(b - a);
    } else {
      a = u;
      u = v;
      v = a + theta2*(b - a);
    }
    cont++;
  }

  // cout << endl << "t: " << (u+v)/2 << endl;
  return (u + v)/2;
}

double Optimization::armijo(double x1, double x2, double d1, double d2, double gama, double eta)
{
  double t = 1;
  double functionStep = function(x1 + t*d1, x2 + t*d2);
  double functionImage = function(x1, x2);
  double regularizer = eta*t*(derivativeX1(x1, x2)*d1 + derivativeX2(x1, x2)*d2);

  // cout << endl << "func step: " << functionStep << endl << "func image: " << functionImage << endl << "reg: " << regularizer << endl; 

  while (functionStep > functionImage + regularizer){
    t = gama*t;
    functionStep = function(x1 + t*d1, x2 + t*d2);
    regularizer = eta*t*(derivativeX1(x1, x2)*d1 + derivativeX2(x1, x2)*d2);

    // cout << endl << "regularizer: " << regularizer << endl;

  }

  return t;
}

double Optimization::gradient(double x1, double x2){
  int k = 0;
  double d1, d2, t, old_x1, old_x2;
  double gradX1 = derivativeX1(x1, x2);
  double gradX2 = derivativeX2(x1, x2);

  double bestX1, bestX2;

  bestX1 = x1;
  bestX2 = x2;
  int it = 0;
  double tolerance = 0.0001;

  while (true){
    d1 = -gradX1;
    d2 = -gradX2;
    t = goldenSection(x1, x2, 0.001, 0.001, d1, d2);//armijo(x1,x2,d1,d2);//0.1;

    old_x1 = x1;
    old_x2 = x2;
    x1 = x1 + t*d1;
    x2 = x2 + t*d2;

    if (x1 < 0 && x2 > 0){
      x1 = x1 + 0.5;
      x2 = x2 - 0.5;
    } else if (x1 > 0 && x2 < 0){
      x1 = x1 - 0.5;
      x2 = x2 + 0.5;
    }
    k++;
    cout << endl << "Gradx1: " << gradX1 << endl << "Gradx2: " << gradX2 << endl << "t: " << t << endl;
    cout << "x1, x2: " << x1 << " " << x2 << endl;

    if (abs(x1 - old_x1) < tolerance && abs(x2 - old_x2) < tolerance) break;

    // if(function(bestX1,bestX2) > function(x1,x2) ){
    //   bestX1 = x1;
    //   bestX2 = x2;
    //   it = k;
    // }

    gradX1 = derivativeX1(x1, x2);
    gradX2 = derivativeX2(x1, x2);
  }

  cout << "Best x: " << x1 << endl << "Best y: " << x2 << endl << "Interacao: " << k << endl;
}

double Optimization::newton(double x1, double x2){
  int k = 0;
  double t, old_x1, old_x2;
  double gradX1 = derivativeX1(x1, x2);
  double gradX2 = derivativeX2(x1, x2);

  double bestX1, bestX2;

  bestX1 = x1;
  bestX2 = x2;
  int it = 0;
  vector<double> d;
  double tolerance = 0.0001;

  cout << "x1: " << x1 << "   x2: " << x2 << endl;
  cout << "gradx1: " << gradX1 << endl;
  cout << "gradx2: " << gradX2 << endl << endl << endl;

  while (true){
    if (k == 10000) break;
    d = dnewton(x1,x2,gradX1,gradX2);
    t = armijo(x1, x2, d[0], d[1]);//goldenSection(x1, x2, 0.001, 0.001, d[0], d[1]);//
    old_x1 = x1;
    old_x2 = x2;
    x1 = x1 + t*d[0];
    x2 = x2 + t*d[1];
    cout << "K: " << k << endl;
    cout << "x1: " << x1 << "   x2: " << x2 << endl;
    cout << "D1: " << d[0] << "   D2: " << d[1] << endl;
    cout << "t: " << t << endl;
    cout << "gradx1: " << gradX1 << endl;
    cout << "gradx2: " << gradX2 << endl << endl << endl;
    k++;
    if (abs(x1) < tolerance && abs(x2) < tolerance){
      cout << "break 1" << endl;
      break;
    }
    if (abs(x1 - old_x1) < tolerance && abs(x2 - old_x2) < tolerance){
      cout << "break 2" << endl;
      break;
    }

    if (x1 < 0 && x2 > 0 && (abs(x1) + abs(x2)) > 2){
      cout << "pai ta aqui" << endl;
      x1 = x1 + 0.7;
      x2 = x2 - 0.7;
    } else if (x1 > 0 && x2 < 0 && (abs(x1) + abs(x2)) > 2){
      cout << "pai ta de celta" << endl;
      x1 = x1 - 0.7;
      x2 = x2 + 0.7;
    }

    // if(function(bestX1,bestX2) > function(x1,x2) ){
    //   bestX1 = x1;
    //   bestX2 = x2;
    //   it = k;
    // }

    gradX1 = derivativeX1(x1, x2);
    gradX2 = derivativeX2(x1, x2);
  }

  cout << endl << "Interacao: " << k << endl;
  cout << "X1: " << x1 << endl << "X2: " << x2 << endl;
};

double Optimization::quasiNewton(double x1, double x2){
  int k = 0;
  double d1, d2, t, old_x1, old_x2;
  vector<vector<double>> H = { { 1, 0}, { 0, 1}};
  vector<vector<double>> invH = { { 1, 0}, { 0, 1}};
  double gradX1 = derivativeX1(x1, x2);
  double gradX2 = derivativeX2(x1, x2);
  double bestX1, bestX2;

  bestX1 = x1;
  bestX2 = x2;
  int it = 0;

  double tolerance = 0.001;

  while (true){
    // if (k == 10) break;
    d1 = - (invH[0][0]*gradX1 + invH[0][1]*gradX2);
    d2 = - (invH[1][0]*gradX1 + invH[1][1]*gradX2);
    t = goldenSection(x1, x2, 0.001, 0.001, d1, d2);//armijo(x1, x2, d1, d2);

    
    old_x1 = x1;
    old_x2 = x2;

    x1 = x1 + t*d1;
    x2 = x2 + t*d2;

    H = BFGS(H,old_x1,old_x2,x1,x2);
    invH = inverse(H);
    k++;
    if (abs(x1) < 0.01 && abs(x2) < 0.2){
      cout << "break" << endl;
      break;
    }
    if (abs(x1 - old_x1) < tolerance && abs(x2 - old_x2) < tolerance){
      cout << "break 2" << endl;
      break;
    }

    if (x1 < 0 && x2 > 0 && (abs(x1) + abs(x2)) > 1.5){
      x1 = x1 + 0.5;
      x2 = x2 - 0.5;
    } else if (x1 > 0 && x2 < 0 && (abs(x1) + abs(x2)) > 1.5){
      x1 = x1 - 0.5;
      x2 = x2 + 0.5;
    }

    gradX1 = derivativeX1(x1, x2);
    gradX2 = derivativeX2(x1, x2);

    cout << "k: " << k << endl;
    cout << "x1: " << x1 << "   x2: " << x2 << endl;
    cout << "D1: " << d1 << "   D2: " << d2 << endl;
    cout << "t: " << t << endl;
    cout << "gradx1: " << gradX1 << "   gradx2: " << gradX2 << endl << endl << endl;
    // cout << H[0][0] << H[0][1] << H[1][0] << H[1][1] << endl;

  }
  cout << "Best x: " << x1 << endl << "Best y: " << x2 << endl << "Interacao: " << k << endl;
}

