#include <iostream>

#include <autodiff/forward/dual.hpp>
using namespace autodiff;

dual f(dual x) {
    return x*x;
}

int main() {
    dual x = 2.0;
    dual u = f(x);

    double dudx = derivative(f, wrt(x), at(x));

    printf("u = %f \n", u);
    printf("du/dx = %f \n", dudx);
}
