#include <stdio.h>
#include <iostream>
#include <math.h> 
#include "otm.h"

using namespace std;

int main()
{
    Optimization opt = Optimization();
    // opt.derivativeX1(1,1);
    // opt.derivativeX2(1,1);
    cout << " Quasi-Newton: " << endl;
    opt.quasiNewton(1, 1);
    return 0;
}