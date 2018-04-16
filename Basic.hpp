#ifndef BASIC_H
#define BASIC_H

#include <ctime>
#include <algorithm>    // for max_element
#include <cmath>        // for sqrt
#include "Matrix.hpp"

namespace Lee{
    
    template<typename T, int N>
    Matrix<T, N, N> eye(){
        Matrix<T, N, N> m;
        for(int i = 0; i != N; ++i)
            m(i, i) = 1;
        return m;
    }

    template<typename T, int N>
    Matrix<T, N, N>& power(Matrix<T, N, N> &m, int k){
        while(--k){
            m = m*m;
        }
        return m;
    }

    template<typename T, int M, int N, int M1, int M2>
    Matrix<T, M, N> row_cat(const Matrix<T, M1, N> &m1, const Matrix<T, M2, N> &m2){
        Matrix<T, M, N> m;
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j){
                if(i < M1) m(i, j) = m1(i, j);
                else m(i, j) = m2(i-M1, j);
            }
        return m;
    }

    template<typename T, int M, int N, int N1, int N2>
    Matrix<T, M, N> col_cat(const Matrix<T, M, N1> &m1, const Matrix<T, M, N2> &m2){
        Matrix<T, M, N> m;
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j){
                if(j < N1) m(i, j) = m1(i, j);
                else m(i, j) = m2(i, j-N1);
            }
        return m;
    }

    template<typename T, int M, int N, int M1>
    Matrix<T, M, N> row_split(const Matrix<T, M1, N> &m, int ri, int rj){
        if(ri>=M1 || rj >= M1) throw std::out_of_range("Matrix index:");

        Matrix<T, M, N> res;
        for(int i = ri; i <= rj; ++i)
            for(int j = 0; j < N; ++j)
                res(i-ri, j) = m(i, j);
        return res;
    }

    template<typename T, int M, int N, int N1>
    Matrix<T, M, N> col_split(const Matrix<T, M, N1> &m, int ci, int cj){
        if(ci>=N1 || cj >= N1) throw std::out_of_range("Matrix index:");

        Matrix<T, M, N> res;
        for(int i = 0; i < M; ++i)
            for(int j = ci; j <= cj; ++j)
                res(i, j-ci) = m(i, j);
        return res;
    }

    template<typename T, int M, int N>
    T max(const Matrix<T, M, N> &m){
        return *std::max_element(m.cbegin(), m.cend());
    }

    template<typename T, int M, int N>
    Matrix<T, M, N> abs(const Matrix<T, M, N> &m){
        Matrix<T, M, N> res = m;
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j)
                res(i, j) = ((m(i, j) > 0) ? m(i, j) : -m(i, j));
        return res;
    }

    template<typename T, int M, int N>
    void Matrix<T, M, N>::permute(int r1, int r2){
        if(r1 >= M || r2 >= M) throw std::out_of_range("Matrix index");
        for(int j = 0; j != N; ++j){
            std::swap((*this)(r1, j), (*this)(r2, j));
        }
    }

    template<typename T, int M, int N>
    Matrix<T, M, N> rand(){
        Matrix<T, M, N> tmp;
        std::srand(time(nullptr));
        for(int i = 0; i != M; ++i)
            for(int j = 0; j != N; ++j)
                tmp(i, j) = std::rand()%20+1;
        return tmp;
    }

    template<typename T, int M1, int N1, int M2, int N2>
    bool is_commute(Matrix<T, M1, N1> &m1, Matrix<T, M2, N2> &m2){
        if(M1 != N2 || M2 != N1) return false;
        if(m1*m2==m2*m1) return true;
        else return false;
    }

    template<typename T, int M, int N>
    int rank(const Matrix<T, M, N> &m){
        return pivot(m).size();
    }

    template<typename T, int M, int N>
    Matrix<T, N, M> transpose(const Matrix<T, M, N> &m){
        Matrix<T, N, M> res;
        for(int i = 0; i != N; ++i)
            for(int j = 0; j != M; ++j)
                res(i, j) = m(j, i);
        return res;
    }

    template<typename T, int M, int N>
    Matrix<T, M-1, N-1> left(const Matrix<T, M, N> &m, int ii, int jj){
        if(M<1 || N<1 || ii<0 || ii>=M || jj<0 || jj>=N) throw std::out_of_range("Matrix index");
        Matrix<T, M-1, N-1> res;
        for(int i = 0; i != M-1; ++i)
            for(int j = 0; j != N-1; ++j){
                if(j>=jj && i<ii) res(i, j) = m(i, j+1);
                else if(i>=ii && j<jj) res(i, j) = m(i+1, j);
                else if(i>=ii && j>=jj) res(i, j) = m(i+1, j+1);
                else res(i, j) = m(i, j);
            }
        return res;
    }

    template<typename T, int N>
    T det(Matrix<T, N, N> m){
        if(N==1) return m(0, 0);
        double num = 0;
        int cntr = 0, cntc = 0;
        for(int i = 0; i != N; ++i) if(!m(i, 0)) ++cntr;
        for(int i = 0; i != N; ++i) if(!m(0, i)) ++cntc;
        if(cntr > cntc)
        {
            for(int i = 0; i != N; ++i){
                if(!m(i, 0)) num += 0;
                else num += (m(i, 0)*cofactor(m, i, 0));
            }
        }
        else{
            for(int i = 0; i != N; ++i){
                if(!m(0, i)) num += 0;
                else num += (m(0, i)*cofactor(m, 0, i));
            }
        }
        return num;
    }

    template<>
    double det<double, 1>(Matrix<double, 1, 1> m){
        return m(0, 0);
    }

    template<typename T, int N>
    T cofactor(Matrix<T, N, N> m, int i, int j){
        return std::pow(-1, i+j)*det(left(m, i, j));
    }

    template<typename T, int N>
    Matrix<T, N, N> adj(const Matrix<T, N, N> &m){
        Matrix<T, N, N> res;
        for(int i = 0; i != N; ++i)
            for(int j = 0; j != N; ++j)
                res(i, j) = cofactor(m, j, i);
        return res;
    }

    template<typename T, int N>
    Matrix<T, N, N> inv(Matrix<T, N, N> &m){
        Matrix<T, N, N> res;        

        if(m.is_invertible()){
            if(m.is_diagonal()){
                for(int i = 0; i < N; ++i)
                    res(i, i) = 1/m(i, i);
            }
            else{
                Matrix<T, N, 2*N> aug = col_cat<T, N, 2*N>(m, eye<T, N>());
                aug = rref(aug);
                res = col_split<T, N, N>(aug, N, 2*N-1);
            }
        }
        return res;
    }

    template<typename T, int N>
    T norm2(const Matrix<T, N, 1> &m){
        return sqrt((transpose(m)*m)(0, 0));
    }
}

#endif