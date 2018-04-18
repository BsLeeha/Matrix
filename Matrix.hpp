#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <array>
#include <cstdlib>
#include <tuple>

// since array provides move constructor and move assignment, we don't need to provide anymore.
// std array doesn't provide any arithmetical operation
// std array do nothing to the default constructor
// We don't need to provide size match test for =, +=, -= for MdArray &a says a is of M*N
// Difficult error to detect: const to non-const
// when writing a reference parameter, think of the rvalue
// ^ is used for bit-wise operation, don't overload it, use power
// we can assign rvalue to *this: (transpose(A+B)==(transpose(A)+transpose(B)))
// nested loop: use < not !=

/**
 * Basic matrix operations: sum, scalar multiple, matrix multiple, commutativity judge(A*B==B*A), transpose, left(i, j)
 * For square matrices: power, determinant, cofactor, adjoint, inverse
**/

// needs to do: special class: vector(N*1), Identity matrix(N*N)

namespace Lee{
    using std::cin;
    using std::cout;
    
    template<typename T, int M, int N>
    class Matrix{
    public:
        // typedef std::pair<int, int> _pos;
        // typedef std::pair<T, _pos> piv;
        // std::vector<piv> pivots;

        typedef std::tuple<int, int, int> piv;
        std::vector<piv> pivots;

        Matrix(){for(auto &i : _data) i = 0;}
        Matrix(std::initializer_list<T> il);
        Matrix(const Matrix &m);              // copy constructor
//        Matrix(Matrix &&m);                 // move constructor
        ~Matrix(){ _data.~array(); }

        T operator()(int i, int j) const { return _data[i*N+j]; }   // data visit: must be const
        T& operator()(int i, int j) { return _data[i*N+j]; }        // for operator >> to write
        Matrix& operator=(const Matrix &m);   // copy assignment
//        Matrix& operator=(Matrix &&m);      // move constructor
        Matrix& operator+=(const Matrix &m);
        Matrix& operator-=(const Matrix &m);
        Matrix& operator*=(double r);
        bool operator==(const Matrix &m) const;
        bool operator!=(const Matrix &m) const;        

        void permute(int r1, int r2);
        bool is_invertible(){ return (M==N) && (rank(*this) == N); }

        void to_eye();
        void to_one();
        void to_zero();

        bool is_diagonal();

         T* begin() { return &_data[0]; }
         T* end() { return &_data[0]+M*N; }

         const T* cbegin() const { return &_data[0]; }
         const T* cend() const {return &_data[0]+M*N; }

        int row() const { return M; }
        int col() const { return N; }

        Matrix<T, M, 1> getcol(int c);
        void setcol(int c, const Matrix<T, M, 1> &);

    private:
        std::array<T, M*N> _data;

    };

    template<typename T, int M, int N>
    std::ostream& operator<<(std::ostream &os, const Matrix<T, M, N> &m){
        for(int i = 0; i < M; ++i){
            for(int j = 0; j < N; ++j){
                os.width(8);
                os.precision(3);
                os << std::fixed << m(i, j) << " ";
            }
            os << "\n";
        }
        return os;
    }

    template<typename T, int M, int N>
    std::istream& operator>>(std::istream &is, Matrix<T, M, N> &m){
        char c;
        if(cin>>c && c=='['){            // start with a '['
            for(int i = 0; i < M; ++i){
                for(int j = 0; j < N; ++j){
                    is >> m(i, j);
                }
                if(i!=M-1 && (is.get()!= ';')) { std::cout <<"\nmissing delimeter ';'\n"; return is; }
            }
            if(is.get()!=']') std::cout << "\nmissing ']', but that's ok!\n";
            return is;
        }
        else std::cout << "\nshould start with '['\n";
        return is;
    }

    template<typename T, int M, int N>
    Matrix<T, M, N>::Matrix(const Matrix &m)     // copy constructor
        : _data{m._data}
    {}

    template<typename T, int M, int N>
    Matrix<T, M, N>::Matrix(std::initializer_list<T> il){
        std::copy(il.begin(), il.end(), _data.begin());
    }

    template<typename T, int M, int N>
    Matrix<T, M, N>& Matrix<T, M, N>::operator=(const Matrix& m){    // copy assignment
//        if(row()!=m.row() || col()!=m.col())
//            throw std::runtime_error("array size mismatch for operator =");
        _data = m._data;
        return *this;
    }

    template<typename T, int M, int N>
    bool Matrix<T, M, N>::operator==(const Matrix &m) const{
        for(int i = 0; i != M; ++i)
            for(int j = 0; j != N; ++j)
                if(m(i, j) != (*this)(i, j)) { cout << "i: " << i <<  " j: " << j << " "<< (m(i, j) == (*this)(i, j)) << "\n";return false;}
        return true;
    }

    template<typename T, int M, int N>
    bool Matrix<T, M, N>::operator!=(const Matrix &m) const{
        return !(*this == m);
    }

    template<typename T, int M, int N>
    Matrix<T, M, N>& Matrix<T, M, N>::operator+=(const Matrix &m){
//        if(row()!=m.row() || col()!=m.col())
//            throw std::runtime_error("array size mismatch for operator +=");
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j)
                (*this)(i, j) += m(i, j);
        return *this;
    }

    template<typename T, int M, int N>
    Matrix<T, M, N>& Matrix<T, M, N>::operator-=(const Matrix &m){
//        if(row()!=m.row() || col()!=m.col())
//            throw std::runtime_error("array size mismatch for operator -=");
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j)
                (*this)(i, j) -= m(i, j);
        return *this;
    }

    template<typename T, int M, int N>
    Matrix<T, M, N>& Matrix<T, M, N>::operator*=(double r){
        for(int i = 0; i != M; ++i)
            for(int j = 0; j != N; ++j)
                (*this)(i, j) *= r;
        return *this;
    }

    template<typename T, int M, int N>
    Matrix<T, M, N> operator+(Matrix<T, M, N> m1, Matrix<T, M, N> m2){
        return m1 += m2;
    }

    template<typename T, int M, int N>
    Matrix<T, M, N> operator-(Matrix<T, M, N> m1, Matrix<T, M, N> m2){
        return m1 -= m2;
    }

   template<typename T, int M, int N>
    Matrix<T, M, N> operator*(double r, Matrix<T, M, N> m){
        return m *= r;
    }

   template<typename T, int M, int N>
    Matrix<T, M, N> operator*(Matrix<T, M, N> m, double r){
        return m *= r;
    }

    template<typename T, int M, int N>
    Matrix<T, M, N> operator/(Matrix<T, M, N> m, double r){
        return m *= (1/r);
    }

    template<typename T, int M1, int N1, int M2, int N2>
    Matrix<T, M1, N2> operator*(const Matrix<T, M1, N1> &m1, const Matrix<T, M2, N2> &m2){
        if(N1 != M2) throw std::runtime_error("array size mismatch for operator *");
        Matrix<T, M1, N2> result;
        for(int i = 0; i < M1; ++i)
            for(int j = 0; j < N2; ++j)
                for(int k = 0; k < N1; ++k)
                  result(i, j) += m1(i, k)*m2(k, j);
        return result;
    }

    template<typename T, int M, int N>
    void Matrix<T, M, N>::to_eye(){
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j){
                if(i == j)(*this)(i, j) = 1;
                else (*this)(i, j) = 0;
            }
    }

    template<typename T, int M, int N>
    void Matrix<T, M, N>::to_one(){
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j)
                (*this)(i, j) = 1;
    }

    template<typename T, int M, int N>
    void Matrix<T, M, N>::to_zero(){
        for(int i = 0; i < M; ++i)
            for(int j = 0; j < N; ++j)
                (*this)(i, j) = 0;
    }

    template<typename T, int M, int N>
    bool Matrix<T, M, N>::is_diagonal(){
        if(M != N) return false;
        int flag = 0;
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < N; ++j)
                if(i!=j && (*this)(i, j)) { flag = 1; break; }
        return flag ? false : true;
    }

    template<typename T, int M, int N>
    Matrix<T, M, 1> Matrix<T, M, N>::getcol(int c){
        Matrix<T, M, 1> col;
        for(int i = 0; i < M; ++i)
            col(i, 0) = (*this)(i, c);
        return col;
    }

    template<typename T, int M, int N, int M1>
    Matrix<T, M1, 1> getcol(Matrix<T, M, N> m, int c, int rs, int re){
        Matrix<T, M1, 1> col;
        if(M1-1 == re-rs) {
            for(int i = rs; i <= re; ++i)
                col(i-rs, 0) = m(i, c);
        }
        return col;
    }

    template<typename T, int M, int N>
    void Matrix<T, M, N>::setcol(int c, const Matrix<T, M, 1> &obj){
        for(int i = 0; i < M; ++i)
            (*this)(i, c) = obj(i, 0);
    }

    template<typename T, int M, int N, int M1>
    void setcol(Matrix<T, M, N> &m, int c, const Matrix<T, M1, 1> &obj, int rs, int re){
        if(M1-1 == re-rs){
            for(int i = rs; i <= re; ++i)
                m(i, c) = obj(i-rs, 0);
        }
    }

    template<typename T, int M, int N>
    Matrix<T, N, 1> least_square(const Matrix<T, M, N> &A, const Matrix<T, M, 1> &b){
        Lee::Matrix<double, N, N> S = transpose(A)*A;
        Lee::Matrix<double, N, 1> x = inv(S)*transpose(A)*b;
        Lee::Matrix<double, M, 1> p = A*x;
        Lee::Matrix<double, M, 1> e = b-p;
        Lee::Matrix<double, M, M> P = A*inv(S)*transpose(A);

        // cout << "A and b \n";
        // cout << A << "\n";
        // cout << b << "\n";
        // cout << "solution to x: \n" << x << "\n";
        // cout << "the error: \n" << e << "\n";
        // cout << "the projection matrix:\n" << P << "\n";
        // cout << "the transpose matrix:\n" << transpose(P) << "\n";
        // cout << std::boolalpha << (P == transpose(P)) << (P*P == P) << (P(0, 0) == transpose(P)(0, 0));
        return x;
    }
}

#endif