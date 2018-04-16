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
}

#endif