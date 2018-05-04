#ifndef _EQUATION_SOLVING_H
#define _EQUATION_SOLVING_H

#include "Polynomial.hpp"

namespace Lee{
    double bisect(const Poly &p, double a, double b, double tol){
        double c;
        if(p(a) == 0) return a;
        if(p(b) == 0) return b;
        if(p(a)*p(b) >= 0) { std::cerr << "f(a)f(b)<0 not satisfied!\n"; return 0; }

        while((b-a)/2 > tol){
            c = (a+b)/2;
            if(p(c) == 0) break;
            if(p(a)*p(c) < 0) b = c;
            else a = c;
        }
        return (a+b)/2;
    }    
}

#endif