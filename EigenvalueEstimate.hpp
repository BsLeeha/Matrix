#ifndef _EIGENVALUEESTIMATE
#define _EIGENVALUEESTIMATE

#include <tuple>
#include "Basic.hpp"

namespace Lee{
    template<typename T, int M, int N>
    std::tuple<T, Matrix<T, M, 1>> power_method(const Matrix<T, M, N> &A, const Matrix<T, M, 1> &x0){
        Matrix<T, M, 1> xp = x0;        // pervious vector
        Matrix<T, M, 1> y = A*xp;
        T u = Lee::max(y);            // eigenvalue                
        Matrix<T, M, 1> xn = y/u;      // next vector
        std::tuple<T, Matrix<T, M, 1>> res;
        double tol = 0.0002;
        int k = 0;

        while(Lee::max(Lee::abs(xn-xp)) > tol){
            std::cout << "x" << k << "\n" << xp;
            std::cout << "Ax" << k << "\n" << xn;
            std::cout << "u" << k << "\n" << "   " << u << "\n\n";            
            xp = xn;
            y = A*xp;
            u = Lee::max(y);
            xn = y/u;
            ++k;
        }
        
        std::get<0>(res) = u;           // eigenvalue
        std::get<1>(res) = xn;          // eigenvector

        std::cout << "largest eigenvalue: " << u << "\n";
        std::cout << "its eigenvector: \n" << xn << "\n";
        return res;
    }
}
#endif

