#ifndef ELIMINATION_H
#define ELIMINATION_H

#include <vector>
#include <tuple>
#include "Matrix.hpp"
#include "Basic.hpp"

namespace Lee{

    // Row echelon form: A -> U
    // Algorithm: Gaussian Elimination
    template<typename T, int M, int N>
    Matrix<T, M, N> upper(Matrix<T, M, N> m){
        
        int flag = 0;
        for (int i = 0; i < M; ++i){        // row pos of pivot
           for (int j = i; j < N; ++j){     // col pos of pivot
                flag = 0;
                if (!m(i, j)){               // pivot zero then exchange rows, no row for exchanging then this col has no pivot
                    for (int r = i+1; r < M; ++r){
                        if (m(r, j)) {m.permute(i, r); flag = 1; break;}  
                    }
                }
                if (m(i, j) || flag){        // pivot not zero(may have some permutations)
                    for (int r = i+1; r < M; ++r){    // forward elimination
                        if (!m(r, j)) continue;
                        T base = m(r, j)/m(i, j);
                        for (int c = j; c < N; ++c){
                            m(r, c) -= (base*m(i, c));
                        }
                    }
                    break;
                }
           }
        }

        return m;
    }

    template<typename T, int M, int N>
    Matrix<T, M, N> lower(Matrix<T, M, N> m){
        int flag = 0;
        for(int i = M-1; i >= 0; --i) {             // row pos of pivot
            for(int j = N-1; j >= 0; --j){          // col pos of pivot
                flag = 0;
                if(!m(i, j)){                 // pivot zero, do permutation
                    for(int r = i-1; r >= 0; --r){
                        if(m(r, j)) { m.permute(i, r); flag = 1; }
                    }
                }
                if(m(i, j) || flag){          // pivot ok
                    for(int r = i-1; r >= 0; --r){      // forward elimination
                        if(!m(r, j)) continue;
                        T base = m(r, j)/m(i, j);
                        for(int c = N-1; c >= 0; --c){
                            m(r, c) -= base*m(i, c);
                        }
                    }
                    break;
                }
            }
        }
        return m;
    }


    template<typename T, int M, int N>
    Matrix<T, M, N> rref(Matrix<T, M, N> m){

        int flag = 0;
        for (int i = 0; i < M; ++i){        // row pos of pivot
           for (int j = i; j < N; ++j){     // col pos of pivot
                flag = 0;
                if (!m(i, j)){               // pivot zero
                    for (int r = i+1; r < M; ++r){
                        if (m(r, j)) {m.permute(i, r); flag = 1; break;}   // exchange rows
                    }
                }
                if (m(i, j) || flag){        // pivot not zero(may have some permutations)
                    for(int c = j; c < N; ++c){       // pivot row turn to identity
                        m(i, c) /= m(i, j);
                    }

                    for (int r = i+1; r < M; ++r){    // forward elimination
                        if (!m(r, j)) continue;
                        T base = m(r, j);
                        for (int c = j; c < N; ++c){
                            m(r, c) -= (base*m(i, c));
                        }
                    }

                    break;
                }
           }
        }

        // for (auto c = m.pivots.rbegin(); c != m.pivots.rend()-1; ++c){      // do back elimination
        //     for(int i = std::get<1>(*c)-1; i >= 0; --i){
        //         T base = m1(i, std::get<2>(*c));       // the elem one row before the pivot
        //         for(int j = 0; j < N; ++j){
        //             m1(i, j) -= (base*m1(std::get<1>(*c), j));
        //         }
        //     }
        // }

        return m;
    }    

    // template<typename T, int M, int N>
    // void solve(Matrix<T, M, N> &A, const Matrix<T, M, 1> &b){
    //     Matrix<T, M, N+1> aug = col_cat<T, M, N+1>(A, b);
    //     aug = rref(aug);
    //     Matrix<T, M, 1> nb = col_split<T, M, 1>(aug, N, N);
    //     Matrix<T, N, 1> res;

    //     cout << "rref(Ab):\n";
    //     cout << aug << "\n";

    //     int zero = 0, one = 0, infty = 0;

    //     if(A.pivots.empty()) rref(A);

    //     for(int i = A.rank(); i != M; ++i)                  // form 0 = b, no solution
    //         if(aug(i, N) != 0) { zero = 1; break; }
    //     if(A.rank()==N && (A.rank() == M || !zero)) one = 1;    
    //     if(!one && !zero) infty = 1;

    //     if(zero){               // for least square, A must have independent cols, or A^TA cannot inverse.
    //         if(A.rank() != N) cout << "Sorry, Least square also cannot solve it!\n";
    //         else{
    //             cout << "Ax=b has no solutions. Here's its best solution:\n";
    //             res = least_square(A, b);
    //             cout << res;
    //         }
    //     }

    //     if(one){
    //         for(int i = 0; i < M; ++i)
    //             res(i, 0) = nb(i, 0);    
    //         if(A.rank()==M)cout << "A is invertible and has one solution: \n";
    //         else cout << "A has " << (M-A.rank()) << " useless rows and has one solution: \n";
    //         cout << res;            
    //     }

    //     if(infty){
    //         cout << "A has " << (N-A.rank()) << " free variables.\n";

    //         if(b != Matrix<T, M, 1>{})
    //         {
    //             cout << "Here's its particular solution: \n";
    //             for(auto c : A.pivots) res(std::get<2>(c) , 0) = nb(std::get<1>(c), 0);
    //             cout << res;
    //         }
    //     }
    // }    
}

#endif