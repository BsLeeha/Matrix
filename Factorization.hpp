#ifndef FACTORIZATION_H
#define FACTORIZATION_H

#include <tuple>
#include "Matrix.hpp"
#include "Basic.hpp"

namespace Lee{
    // If Permutation is done on A, then E is product of Es and P, so L may not be lower triangular.
    template<typename T, int N> 
    std::tuple<Matrix<T, N, N>, Matrix<T, N, N>> PLU(Matrix<T, N, N> A){
        Matrix<T, N, N> L;
        Matrix<T, N, N> E, tmp;
        Matrix<T, N, N> P;
        E.to_eye(); tmp.to_eye(); P.to_eye(); L.to_eye();
        std::tuple<Matrix<T, N, N>, Matrix<T, N, N>> res;

        for(int i = 0; i < N; ++i)          // row pos of pivots
            for(int j = i; j < N; ++j){     // col pos of pivots 
                if(!A(i, j)) {                         // bad pivot, do permutations
                    for(int r = i+1; r < N; ++r){
                        if(A(r, j)) { 
                            P.to_eye();
                            A.permute(r, i); 
                            P(r, r) = 0;
                            P(i, i) = 0;
                            P(r, i ) = 1;
                            P(i, r) = 1; 
                            break;
                        }
                    }
                }
                if(A(i, j)){                    // good pivot, do forward elimination
                    tmp.to_eye();
                    for(int r = i+1; r < N; ++r){
                        if(!A(r, j)) continue;
                        T base = A(r, j)/A(i, j);           // keng!!!
                        tmp(r, j) =  -base;             
                        for(int c = j; c < N; ++c){
                            A(r, c) -= (base*A(i, c));
                        }
                    }
                    E = tmp*P*E;
                    break;                              // find pivot in the next line
                }
                else continue;
            }

        L = inv(E);
        std::get<0>(res) = L;
        std::get<1>(res) = A;
        return res;
    }

    template<typename T, int M, int N>
    std::tuple<Matrix<T, M, N>, Matrix<T, N, N>> QRGramScmidt(Matrix<T, M, N> A){
        Matrix<T, M, N> Q;
        Matrix<T, N, N> R;
        std::tuple<Matrix<T, M, N>, Matrix<T, N, N>> res; 
        Matrix<T, M, 1> y;
        double tmp;

        for(int i = 0; i < N; ++i){
            y.to_zero();
            y = A.getcol(i);
            for(int j = 0; j < i; ++j){
                tmp = (transpose(Q.getcol(j))*A.getcol(i))(0, 0);
                y -= Q.getcol(j)*(tmp);
                R(j, i) = tmp;
            }
            R(i, i) = norm2(y);
            Q.setcol(i, y/norm2(y));
        }

        std::get<0>(res) = Q;
        std::get<1>(res) = R;
        return res;
    }

    template<typename T, int M, int N>
    std::tuple<Matrix<T, M, M>, Matrix<T, M, N>> QRGramScmidtex(Matrix<T, M, N> A){
        Matrix<T, M, M> Q;
        Matrix<T, M, N> R;
        std::tuple<Matrix<T, M, M>, Matrix<T, M, N>> res; 
        Matrix<T, M, 1> y;
        double tmp;

        Matrix<T, M, M> Aex;
        for(int i = 0; i < N; ++i)
            Aex.setcol(i, A.getcol(i));
        do{
            for(int i = N; i < M; ++i)
                Aex.setcol(i, rand<T, M, 1>());
        }while(rank(Aex) != M);

        for(int i = 0; i < M; ++i){
            y = Aex.getcol(i);
            for(int j = 0; j < i; ++j){
                tmp = (transpose(Q.getcol(j))*Aex.getcol(i))(0, 0);
                y -= Q.getcol(j)*tmp;
                if(i < N)R(j, i) = tmp;
            }
            R(i, i) = norm2(y);
            Q.setcol(i, y/norm2(y));
        }

        std::get<0>(res) = Q;
        std::get<1>(res) = R;
        return res;
    }    

}

#endif