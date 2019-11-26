#include <stdio.h>
#include <iostream>
#include <math.h> 
#include "otm.h"

int main()
{
    Optimization opt = Optimization();
    // opt.derivativeX1(1,1);
    // opt.derivativeX2(1,1);
    opt.gradient(1, 1);
    return 0;
}