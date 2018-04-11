#ifndef FACTORIZATION_H
#define FACTORIZATION_H

#include <tuple>
#include "Matrix.hpp"
#include "Basic.hpp"

namespace Lee{
    // If Permutation is done on A, then E is product of Es and P, so L may not be lower triangular.
    template<typename T, int N> 
    std::tuple<Matrix<T, N, N>, Matrix<T, N, N>> PLU(const Matrix<T, N, N> &m){
        Matrix<T, N, N> L;
        Matrix<T, N, N> U = m;
        Matrix<T, N, N> E, tmp;
        Matrix<T, N, N> P;
        E.to_eye(); tmp.to_eye(); P.to_eye();
        std::tuple<Matrix<T, N, N>, Matrix<T, N, N>> res;

        int flag = 0;
        for(int i = 0; i < N; ++i)          // row pos of pivots
            for(int j = i; j < N; ++j){     // col pos of pivots 
                flag = 0;
                if(!U(i, j)) {                         // bad pivot, do permutations
                    for(int r = i+1; r < N; ++r){
                        if(U(r, j)) { 
                            P.to_eye();
                            U.permute(r, i); 
                            flag = 1;
                            P(r, r) = 0;
                            P(i, i) = 0;
                            P(r, i ) = 1;
                            P(i, r) = 1; 
                            break;
                        }
                    }
                }
                if(U(i, j) || flag){                    // good pivot, do forward elimination
                    tmp.to_eye();
                    for(int r = i+1; r < N; ++r){
                        if(!U(r, j)) continue;
                        T base = U(r, j)/U(i, j);           // keng!!!
                        tmp(r, j) =  -base;             
                        for(int c = j; c < N; ++c){
                            U(r, c) -= (base*U(i, c));
                        }
                    }

                    E = tmp*P*E;
                    break;                              // find pivot in the next line
                }
            }
        L = inv(E);
        std::get<0>(res) = L;
        std::get<1>(res) = U;
        return res;
    }

}

#endif