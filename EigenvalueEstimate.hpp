#ifndef _EIGENVALUEESTIMATE
#define _EIGENVALUEESTIMATE

#include <tuple>
#include "Basic.hpp"

namespace Lee{
    // template<typename T, int M, int N>
    // std::tuple<T, Matrix<T, M, 1>> power_method(const Matrix<T, M, N> &A, const Matrix<T, M, 1> &x0, T tol){
    //     Matrix<T, M, 1> xp = x0;        // pervious vector
    //     Matrix<T, M, 1> y = A*xp;
    //     T u = Lee::max(y);            // eigenvalue                
    //     Matrix<T, M, 1> xn = y/u;      // next vector
    //     cout << y << "\n" << u;
    //     std::tuple<T, Matrix<T, M, 1>> res;
    //     int k = 0;

    //     cout << xn << "\n" << xp;
    //     while(Lee::max(Lee::abs(xn-xp)) > tol){
    //         std::cout << "x" << k << "\n" << xp;
    //         std::cout << "Ax" << k << "\n" << xn;
    //         std::cout << "u" << k << "\n" << "   " << u << "\n\n";            
    //         xp = xn;
    //         y = A*xp;
    //         u = Lee::max(y);
    //         xn = y/u;
    //         ++k;
    //     }
        
    //     std::get<0>(res) = u;           // eigenvalue
    //     std::get<1>(res) = xn;          // eigenvector

    //     std::cout << "largest eigenvalue: " << u << "\n";
    //     std::cout << "its eigenvector: \n" << xn << "\n"; 
    //     return res;
    // }

    template<typename T, int N>
    std::tuple<T, Matrix<T, N, 1>> power_method(const Matrix<T, N, N> &A, const Matrix<T, N, 1> &x0){
        Matrix<T, N, 1> u;
        Matrix<T, N, 1> x = x0;
        std::tuple<T, Matrix<T, N, 1>> res;
        T lambda;

        for(int i = 0; i < 40; ++i){
            u = x/norm2(x);
            x = A*u;
            lambda = (transpose(u)*x)(0, 0);
        }

        std::get<0>(res) = lambda;
        std::get<1>(res) = x;
        return res;
    }

    template<typename T, int N>
    std::tuple<T, Matrix<T, N, 1>> inverse_power_method(const Matrix<T, N, N> &A, const Matrix<T, N, 1> &x0, double s){
        Matrix<T, N, 1> u;
        Matrix<T, N, 1> x = x0;
        double lambda;
        std::tuple<T, Matrix<T, N, 1>> res;

        for(int i = 0; i < 40; ++i){
            u = x/norm2(x);
            x = GaussianDirect(A-s*eye<T, N>(), u);
            lambda = (transpose(u)*x)(0, 0); 
        }
        lambda = 1/lambda + s;

        std::get<0>(res) = lambda;
        std::get<1>(res) = x;
        return res;
    }
}
#endif

