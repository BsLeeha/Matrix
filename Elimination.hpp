#ifndef ELIMINATION_H
#define ELIMINATION_H

#include <vector>
#include <tuple>
#include "Matrix.hpp"

namespace Lee{

    template<typename T>
    std::vector<T> pivot(Matrix<T, 2> m){
        size_t M = m.rows();
        size_t N = m.cols();

        std::vector<T> pivots;
        T piv;

        for(size_t i = 0; i < m.rows(); ++i)                      // row pos of pivot
            for(size_t j = i; j < m.cols(); ++j){                 // col pos of pivot
                if(!m(i, j)){                           // pivot zero, row permutation
                    for(size_t r = i+1; r < M; ++r){
                        if(m(r, j)) { m.rowPermute(i, r); break; }
                    }
                }
                if(m(i, j)){                    // pivot not zero, forward elimination
                    piv = m(i, j);
                    pivots.push_back(piv);
                    for(size_t r = i+1; r < M; ++r){
                        if(!m(r, j)) continue;          // variable zero, no need to eliminate
                        T base = m(r, j)/m(i, j);
                        for(size_t c = j; c < N; ++c){
                            m(r, c) -= (base*m(i, c));
                        }
                    }
                    break;                              // pivot find in this col, break
                }
                else continue;                          // zero col, find pivot in next col
            }
        return pivots;
    }

    // Row echelon form: A -> U
    // Algorithm: Gaussian Elimination
    template<typename T>
    Matrix<T, 2> upper(Matrix<T, 2> m){
        size_t M = m.rows();
        size_t N = m.cols();

        for (size_t i = 0; i < M; ++i){                 // row pos of pivot
           for (size_t j = i; j < N; ++j){              // col pos of pivot
                if (!m(i, j)){                       // pivot zero, row permutation
                    for (size_t r = i+1; r < M; ++r){
                        if (m(r, j)) {m.rowPermute(i, r); break;}  
                    }
                }
                if (m(i, j)){                // pivot not zero, forward elimination
                    for (size_t r = i+1; r < M; ++r){  
                        if (!m(r, j)) continue;      // variable zero, no need to eliminate
                        T base = m(r, j)/m(i, j);
                        for (size_t c = j; c < N; ++c){
                            m(r, c) -= (base*m(i, c));
                        }
                    }
                    break;                          // pivot find in this col, break
                }
                else continue;                      // zero col, find pivot in next col
           }
        }

        return m;
    }

    template<typename T>
    Matrix<T, 2> lower(Matrix<T, 2> m){
        size_t M = m.rows();
        size_t N = m.cols();

        for(size_t i = M; i-- > 0; ) {             // row pos of pivot
            for(size_t j = N; j-- > 0; ){          // col pos of pivot
                if(!m(i, j)){                       // pivot zero, row permutation
                    for(size_t r = i; r-- > 0; ){
                        if(m(r, j)) { m.rowPermute(i, r); }
                    }
                }
                if(m(i, j)){                // pivot not zero, forward elimination
                    for(size_t r = i; r-- > 0; ){    
                        if(!m(r, j)) continue;      // variable zero, no need to eliminate
                        T base = m(r, j)/m(i, j);
                        for(size_t c = N; c-- > 0; ){
                            m(r, c) -= base*m(i, c);
                        }
                    }
                    std::cout << m;
                    break;                          // pivot find in this col, break
                }
                else continue;                      // zero col, find pivot in next col
            }
        }
        return m;
    }


    template<typename T>
    Matrix<T, 2> rref(Matrix<T, 2> m){
        size_t M = m.rows();
        size_t N = m.cols();

        std::vector<std::tuple<T, int, int>> pivots;
        std::tuple<T, int, int> pivot;        

        for (size_t i = 0; i < M; ++i){                // row pos of pivot
           for (size_t j = i; j < N; ++j){             // col pos of pivot
                if (!m(i, j)){                      // pivot zero, row permutation
                    for (size_t r = i+1; r < M; ++r){
                        if (m(r, j)) {m.rowPermute(i, r); break;}  
                    }
                }
                if (m(i, j)){                       // pivot not zero, forward elimination
                    for(size_t c = j+1; c < N; ++c){     // pivot row turn to identity
                        m(i, c) /= m(i, j);
                    }
                    std::get<0>(pivot) = m(i, j);
                    std::get<1>(pivot) = i;
                    std::get<2>(pivot) = j;
                    pivots.push_back(pivot);
                    m(i, j) = 1;                    
                    for (size_t r = i+1; r < M; ++r){    // forward elimination
                        if (!m(r, j)) continue;
                        T base = m(r, j);
                        for (size_t c = j; c < N; ++c){
                            m(r, c) -= (base*m(i, c));
                        }
                    }
                    break;                          // pivot find in this col, break
                }
                else continue;                       // zero col, find pivot in next col
           }
        }

        for (auto p = pivots.rbegin(); p != pivots.rend()-1; ++p){      // do back elimination
            for(size_t r = std::get<1>(*p); r > 0; --r){
                T base = m(r-1, std::get<2>(*p));         
                for(size_t c = r-1; c < N; ++c){
                    m(r-1, c) -= (base*m(std::get<1>(*p), c));
                }
            }
        }

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