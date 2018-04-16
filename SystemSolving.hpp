#ifndef _SYSTEMSOLVING_H
#define _SYSTEMSOLVING_H

#include "Elimination.hpp"
#include "Factorization.hpp"

namespace Lee{

    template<typename T, int N>
    Matrix<T, N, 1> UpperBackSub(const Matrix<T, N, N> &U, Matrix<T, N, 1> b){
        for(int i = N-1; i >= 0; --i){
            for(int j = N-1; j > i; --j)
                b(i, 0) -= (U(i, j)*b(j, 0));
            b(i, 0) /= U(i, i);
        }
        return b;
    }

    template<typename T, int N>
    Matrix<T, N, 1> LowerBackSub(const Matrix<T, N, N> &L, Matrix<T, N, 1> b){
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < i; ++j)
                b(i, 0) -=(L(i, j)*b(j, 0)); 
            b(i, 0) /= L(i, i);
        }
        return b;
    }
    template<typename T, int N>
    Matrix<T, N, 1> GaussianDirect(const Matrix<T, N, N> &A, const Matrix<T, N, 1> &b){
        // Forward elimination
        Matrix<T, N, N+1> Au = col_cat<T, N, N+1>(A, b);
        Au = upper(Au);
        Matrix<T, N, N> An = col_split<T, N, N>(Au, 0, N-1);
        Matrix<T, N, 1> bn = col_split<T, N, 1>(Au, N, N);
        // Back substitution
        Matrix<T, N, 1> x = UpperBackSub(An, bn);

        return x;
    }

    template<typename T, int N>
    Matrix<T, N, 1> GaussianPLU(Matrix<T, N, N> A, Matrix<T, N, 1> b){
        // Forward elimination
        Matrix<T, N, N> L = std::get<0>(PLU(A));
        Matrix<T, N, N> U = std::get<1>(PLU(A));
        // Back substitution
        Matrix<T, N, 1> y = LowerBackSub(L, b);
        Matrix<T, N, 1> x = UpperBackSub(U, y);

        return x;
    }

    // Jacobi iteration in system form
    template<typename T, int N>
    Matrix<T, N, 1> DirectJacobi(Matrix<T, N, N> A, Matrix<T, N, 1> b, Matrix<T, N, 1> x0, T tol){
        Matrix<T, N, 1> x1 = x0, x2;
        int k = 0, km = 20;

        while(k++<km && norm2(A*x2-b)>tol){
            x2.to_zero();
            for(int i = 0; i < N; ++i){
                for(int j = 0; j < i; ++j)
                    x2(i, 0) += A(i, j)*x1(j, 0);
                for(int j = i+1; j < N; ++j)
                    x2(i, 0) += A(i, j)*x1(j, 0);
                x2(i, 0) = -(x2(i, 0)-b(i, 0))/A(i, i);
            }
            cout << "x:\n" << x1;             
            x1 = x2;
        }
        if(k >= km) std::cerr << "Iteration Fail!\n";
        return x2;
    }

    // Jacobi iteration in matrix form
    template<typename T, int N>
    Matrix<T, N, 1> MatrixJacobi(Matrix<T, N, N> A, Matrix<T, N, 1> b, Matrix<T, N, 1> x0, T tol){
        Matrix<T, N, 1> x1 = x0, x2;
        int k = 0, km = 20;
        Matrix<T, N, N> D, L, U;

        for(int i = 0; i < N; ++i)
                D(i, i) = A(i, i);
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < i; ++j)
                L(i, j) = A(i, j);
        for(int i = 0; i < N; ++i)
            for(int j = i+1; j < N; ++j)
                U(i, j) = A(i, j);

        while(k++<km && norm2(A*x2-b)>tol){
            std::cout << "x:\n" << x1;
            x2 = inv(D)*(b-(L+U)*x1);
            x1 = x2;
        }
        if(k >= km) std::cerr << "Iteration Fail!\n";
        return x2;
    }

    template<typename T, int N>
    Matrix<T, N, 1> DirectGaussSeidel(Matrix<T, N, N> A, Matrix<T, N, 1> b, Matrix<T, N, 1>x0, T tol){
        Matrix<T, N, 1> x1 = x0, x2;
        int k = 0, km = 20;

        while(k++<km && norm2(A*x2-b)>tol){
            x2.to_zero();
            for(int i = 0; i < N; ++i){
                for(int j = 0; j < i; ++j)
                    x2(i, 0) += A(i, j)*x2(j, 0);
                for(int j = i+1; j < N; ++j)
                    x2(i, 0) += A(i, j)*x1(j, 0);
                x2(i, 0) = -(x2(i, 0)-b(i, 0))/A(i, i);
            }
            std::cout << "x:\n" << x1;
            x1 = x2;
        }
        if(k >= km) std::cerr << "Iteration Fail!\n";
        return x2;
    }

    // maybe not exist!!!
    template<typename T, int N>
    Matrix<T, N, 1> MatrixGaussSeidel(Matrix<T, N, N> A, Matrix<T, N, 1> b, Matrix<T, N, 1>x0, T tol){
        Matrix<T, N, 1> x1 = x0, x2;
        int k = 0, km = 20;
        Matrix<T, N, N> L, D, U;

        for(int i = 0; i < N; ++i)
            D(i, i) = A(i, i);
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < i; ++j)
                L(i, j) = A(i, j);
        for(int i = 0; i < N; ++i)
            for(int j = i+1; j < N; ++j)
                U(i, j) = A(i, j);

        while(k++<km && norm2(A*x2-b)>tol){
            x2 = inv(D)*(b-U*x1-L*x2);
            x1 = x2;
        }
        if(k >= km) std::cerr<<"Iteration Fail\n";
        return x2;
    }
}

#endif