#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <type_traits>
#include <functional>

/*
** Operation define:
** classes
** | Matrix<T, M, N>
** | Scalar = Matrix<T, 1, 1>
** | RealScalar = Matrix<double, 1, 1>
** | Vetcor = Matrix<T, M, 1> or Matrix<T, 1, N>
** | RealVetcor = Matrix<double, M, 1> or Matrix<T, 1, N>
** | EyeMatrix = Matrix<T, N, N>
** | RealEyeMatrix = Matrix<double, N, N>
** initialize
** | Matrix a
** | Matrix a{...}
** | Matrix a = {...}
** | Matrix a(b)
** | Matrix a = Matrix b
** visit
** | a(...)
** | a.row(i)
** | a.col(i)
** 
*/

namespace MatrixImpl{
    template<typename T, size_t N>
    struct NestedInitializerList;    

    template<typename M>
    void index_bounds_check(const M &m, size_t r, size_t c);

    template<typename M>
    void matrix_valid(const M &m);

    template<typename M, typename...Dims>
    void matrix_valid(const M &m, const Dims &...dims);      

    template<typename T>
    struct IsParenType;          
}

namespace Lee{

    template<typename T, size_t N>
    using NestedInitializerListN = typename MatrixImpl::NestedInitializerList<T, N>::type;

    template<typename T, typename V, typename F>
    class applyProxy;

    template<typename T, typename V1, typename V2, typename F>
    class binaryProxy;    

    struct slice{
        explicit slice(size_t nstart, size_t nend, size_t nstride) 
            : start{nstart*nstride}, size{nend-nstart+1}, stride{nstride} {}

        size_t operator()(size_t i) const{
            return start+i*stride;
        }

        size_t start;
        size_t size;
        size_t stride;
    };
    template<typename T, size_t M, size_t N>
    class sliceMatrix;

    template<typename T, size_t M, size_t N, typename V = std::vector<T>>
    class Matrix{
    public:
        using value_type     = T;
        using iterator       = typename V::iterator;
        using const_iterator = typename V::const_iterator;

        // Constructors
        Matrix(){
            elems.resize(M*N); 
            std::fill(elems.begin(), elems.end(), 0);
            
            assert(elems.size() == M*N && "construction fail");
        }                        

        explicit Matrix(size_t cnt, T elem){
            assert(cnt == M*N && "overinput");

            elems.resize(M*N);
            std::fill(elems.begin(), elems.end(), elem);

            assert(elems.size() == M*N && "construction fail");
        }

        Matrix(NestedInitializerListN<T, 1> il)
        {
            assert(il.size() <= size() && "overinput");            
            
            elems.assign(il);
            elems.insert(elems.end(), size()-elems.size(), static_cast<T>(0));

            assert(elems.size() == M*N && "construction fail");            
        }

        Matrix(NestedInitializerListN<T, 2> il){
            assert(il.size() <= M && "rows overinput");

            for(auto i : il){
                assert(i.size() <= N && "cols overinput");
                
                std::copy(i.begin(), i.end(), std::back_inserter(elems));
                elems.insert(elems.end(), N-i.size(), static_cast<T>(0));
            }
            elems.insert(elems.end(), size()-elems.size(), static_cast<T>(0));

            assert(elems.size() == M*N && "construction fail");                
        }        

        // for several kinds of expression proxies
        Matrix(const V &vec) : elems{vec} {}

        template<typename V1>
        Matrix(const Matrix<T, M, N, V1> &rhs){
            MatrixImpl::matrix_valid(rhs);

            elems.resize(size());
            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) = rhs(i, j); 

            assert(elems.size() == M*N && "assignment fail");                                          
        }

        template<size_t M1, size_t N1>
        Matrix(const sliceMatrix<T, M1, N1> &rhs){
            assert(rhs.rows() == M && rhs.cols() == N && "assignment does not match");

            elems.resize(size());
            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) = rhs(i, j);

            assert(elems.size() == M*N && "assignment fail");                      
        }

        Matrix(Matrix &&rhs) = default;      

        ~Matrix() = default;
        
        // assignments
        Matrix& operator=(const Matrix &rhs){
            if(this != &rhs) elems = rhs.elems;

            return *this;
        }

        Matrix& operator=(const Matrix &&rhs){           
            if(this != &rhs) elems = std::move(rhs.elems);

            return *this;
        }

        Matrix& operator=(NestedInitializerListN<T, 1> il){
            assert(il.size() <= size() && "overassign");            
            
            elems.assign(il);
            elems.insert(elems.end(), size()-elems.size(), static_cast<T>(0));

            assert(elems.size() == M*N && "assignment fail");                       
            return *this;         
        }

        Matrix& operator=(NestedInitializerListN<T, 2> il){
            assert(il.size() <= M && "rows overassign");

            elems.clear();
            
            for(auto i : il){
                assert(i.size() <= N && "cols overassign");
                
                std::copy(i.begin(), i.end(), std::back_inserter(elems));
                elems.insert(elems.end(), N-i.size(), static_cast<T>(0));
            }

            elems.insert(elems.end(), size()-elems.size(), static_cast<T>(0));

            assert(elems.size() == M*N && "assignment fail");  
            return *this;                              
        }        

        template<typename V1>
        Matrix& operator=(const Matrix<T, M, N, V1> &rhs){
            MatrixImpl::matrix_valid(rhs);

            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) = rhs(i, j);
                    
            assert(elems.size() == M*N && "assignment fail");                      
            return *this;
        }

        template<size_t M1, size_t N1>
        Matrix& operator=(const sliceMatrix<T, M1, N1> &rhs){
            assert(rhs.rows() == M && rhs.cols() == N && "assignment does not match");

            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) = rhs(i, j);

            assert(elems.size() == M*N && "assignment fail");                      
            return *this;
        }

        // member functions
        static constexpr size_t rows() { return M; }

        static constexpr size_t cols() { return N; }

        static constexpr size_t size() { return M*N; }

        template<typename F>
        Matrix<T, M, N, applyProxy<T, V, F>> apply(F f){
            return applyProxy<T, V, F>(data(), f);
        }
        
        // element access
        template<typename Q = V>
        typename std::enable_if<MatrixImpl::IsParenType<Q>::value, T&>::type
        operator()(size_t i, size_t j){
            MatrixImpl::index_bounds_check(*this, i, j);   
            
            return elems(i, j);
        }

        template<typename Q = V>
        typename std::enable_if<!MatrixImpl::IsParenType<Q>::value, T&>::type
        operator()(size_t i, size_t j) { 
            MatrixImpl::index_bounds_check(*this, i, j);   

            return elems[i*N+j]; 
        }    

        template<typename Q = V>
        typename std::enable_if<MatrixImpl::IsParenType<Q>::value, T>::type
        operator()(size_t i, size_t j) const {
            MatrixImpl::index_bounds_check(*this, i, j);   
            
            return elems(i, j);
        }

        template<typename Q = V>
        typename std::enable_if<!MatrixImpl::IsParenType<Q>::value, T>::type
        operator()(size_t i, size_t j) const { 
            MatrixImpl::index_bounds_check(*this, i, j);

            return elems[i*N+j]; 
        }  

        V& data() { return elems; }

        const V& data() const { return elems; }

        iterator begin() { return elems.begin(); }

        const_iterator begin() const { return elems.begin(); }

        iterator end() { return elems.end(); }

        const_iterator end() const { return elems.end(); }

        // submatrix access
        sliceMatrix<T, M, N> operator()(slice r, slice c){
            return sliceMatrix<T, M, N>(*this, r, c);
        }
        
        sliceMatrix<T, M, N> row(size_t i){
            MatrixImpl::index_bounds_check(*this, i, 0);

            return sliceMatrix<T, M, N>(*this, slice(i, i, N), slice(0, N-1, 1));
        }

        sliceMatrix<T, M, N> col(size_t i){
            MatrixImpl::index_bounds_check(*this, 0, i);

            return sliceMatrix<T, M, N>(*this, slice(0, M-1, N), slice(i, i, 1));
        }

        // unary operations
        Matrix<T, M, N, applyProxy<T, V, std::negate<T>>> operator-(){
            return apply(std::negate<T>());
        }

        // binary operations
        template<typename V1>
        Matrix<T, M, N, binaryProxy<T, V, V1, std::plus<T>>> operator+(const Matrix<T, M, N, V1> &rhs){
            return binaryProxy<T, V, V1, std::plus<T>>(elems, rhs.data(), std::plus<T>());
        }

        template<typename V1>
        Matrix<T, M, N, binaryProxy<T, V, V1, std::minus<T>>> operator-(const Matrix<T, M, N, V1> &rhs){
            return binaryProxy<T, V, V1, std::minus<T>>(elems, rhs.data(), std::minus<T>());
        }

        // arithmetic operations
        Matrix& operator+=(const Matrix &m){
            MatrixImpl::matrix_valid(*this, m);       
            
            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) += m(i, j);
            return *this;            
        }

        Matrix& operator-=(const Matrix &m){
            MatrixImpl::matrix_valid(*this, m);       
            
            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) -= m(i, j);
            return *this;            
        }

        Matrix& operator*=(double r){
            MatrixImpl::matrix_valid(*this);       
            
            for(size_t i = 0; i != M; ++i)
                for(size_t j = 0; j != N; ++j)
                    (*this)(i, j) *= r;
            return *this;            
        }

        Matrix& operator/=(double r){
            MatrixImpl::matrix_valid(*this);       
            
            return (*this) *= 1/r;
        }

        void to_eye(){
            MatrixImpl::matrix_valid(*this);       
            
            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j){
                    if(i == j)(*this)(i, j) = 1;
                    else (*this)(i, j) = 0;
                }            
        }

        void to_one(){
            MatrixImpl::matrix_valid(*this);       

            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) = 1;            
        }

        void to_zero(){
            MatrixImpl::matrix_valid(*this);       
            
            for(size_t i = 0; i < M; ++i)
                for(size_t j = 0; j < N; ++j)
                    (*this)(i, j) = 0;            
        }

        bool is_diagonal(){
            MatrixImpl::matrix_valid(*this);       
            
            if(M != N) return false;
            int flag = 0;
            for(size_t i = 0; i < N; ++i)
                for(size_t j = 0; j < N; ++j)
                    if(i!=j && (*this)(i, j)) { flag = 1; break; }
            return flag ? false : true;            
        }

    private:
        V elems;

    };

    template<typename T, size_t M, size_t N, typename V>
    std::ostream& operator<<(std::ostream &os, const Matrix<T, M, N, V> &m){
        MatrixImpl::matrix_valid(m);       

        for(size_t i = 0; i < M; ++i){
            for(size_t j = 0; j < N; ++j){
                os.width(8);
                os.precision(3);
                os << std::fixed << m(i, j) << " ";
            }
            os << '\n';
        }
        os << '\n';
        return os;
    }

    template<typename T, size_t M, size_t N>
    std::istream& operator>>(std::istream &is, Matrix<T, M, N> &m){
        char c;
        if(is>>c && c=='['){            // start with a '['
            for(size_t i = 0; i < M; ++i){
                for(size_t j = 0; j < N; ++j){
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

    // concatenate 
    // vertcat
    template<typename T, size_t M1, size_t M2, size_t N>
    Matrix<T, M1+M2, N> vertcat(const Matrix<T, M1, N> &m1, const Matrix<T, M2, N> &m2){
        MatrixImpl::matrix_valid(m1, m2);

        Matrix<T, M1+M2, N> res;
        for(size_t i = 0; i < M1+M2; ++i)
            for(size_t j = 0; j < N; ++j)
                if(i < M1) res(i, j) = m1(i, j);
                else res(i, j) = m2(i-M1, j);

        return res;
    }

    template<typename T, size_t M, size_t N1, size_t N2>
    Matrix<T, M, N1+N2> horzcat(const Matrix<T, M, N1> &m1, const Matrix<T, M, N2> &m2){
        MatrixImpl::matrix_valid(m1, m2);

        Matrix<T, M, N1+N2> res;
        for(size_t i = 0; i < M; ++i)
            for(size_t j = 0; j < N1+N2; ++j)
                if(j < N1) res(i, j) = m1(i, j);
                else res(i, j) = m2(i, j-N1);

        return res;
    }

    template<typename T, typename P>
    struct Iterator{
        const P &proxy;
        size_t pos;

        Iterator(const P &nproxy, size_t npos) : proxy{nproxy}, pos{npos} {}

        const T operator*() {
            return proxy[pos];
        }

        Iterator& operator++(){
            ++pos;
            return *this;
        }

        Iterator& operator--(){
            --pos;
            return *this;
        }

        Iterator& operator+(size_t i){
            pos += i;
            return *this;
        }

        Iterator& operator-(size_t i){
            pos -= i;
            return *this;
        }

        bool operator==(const Iterator &rhs) const{
            return pos == rhs.pos;
        }

        bool operator!=(const Iterator &rhs) const{
            return pos == rhs.pos;
        }

    };


    template<typename T, typename P>
    struct ConstIterator{
        const P &proxy;
        size_t pos;

        ConstIterator(const P &nproxy, size_t npos) : proxy{nproxy}, pos{npos} {}

        const T operator*() const{
            return proxy[pos];
        }

        ConstIterator& operator++(){
            ++pos;
            return *this;
        }

        ConstIterator& operator--(){
            --pos;
            return *this;
        }

        ConstIterator& operator+(size_t i){
            pos += i;
            return *this;
        }

        ConstIterator& operator-(size_t i){
            pos -= i;
            return *this;
        }

        bool operator==(const ConstIterator &rhs) const{
            return pos == rhs.pos;
        }

        bool operator!=(const ConstIterator &rhs) const{
            return pos == rhs.pos;
        }

    };    

    template<typename T>
    class Scalar{
    public:
        explicit Scalar(const T &e) : elem{e} {}

        const T& operator[](size_t i) const{
            return elem;
        }

    private: 
        T elem;
    };

    template<typename T, typename V, typename F>
    class applyProxy{
    public:
        using value_type     = T;
        using iterator       = Iterator<T, applyProxy>;
        using const_iterator = ConstIterator<T, applyProxy>;

        applyProxy(const V &left, F &function) : lhs{left}, func{function} {}

        T operator[](size_t i) const{
            return func(lhs[i]);
        }

        size_t size() const{
            return lhs.size();
        }

        iterator begin(){
            return Iterator<T, applyProxy>(*this, 0);
        }

        const_iterator begin() const {
            return ConstIterator<T, applyProxy>(*this, 0);            
        }

        iterator end(){
            return Iterator<T, applyProxy>(*this, size());
        }

        const_iterator end() const {
            return ConstIterator<T, applyProxy>(*this, size());
        }

    private:
        const V &lhs;
        const F &func;
    };

    template<typename T, typename V1, typename V2, typename F>
    class binaryProxy{
    public:
        using value_type     = T;
        using iterator       = Iterator<T, binaryProxy>;
        using const_iterator = ConstIterator<T, binaryProxy>;

        binaryProxy(const V1 &left, const V2 &right, const F &f) 
            : lhs{left}, rhs{right}, func{f} {}
        
        T operator[](size_t i) const{
            return func(lhs[i], rhs[i]);
        }

        size_t size() const { return lhs.size(); }

        iterator begin(){
            return Iterator<T, binaryProxy>(*this, 0);
        }

        const_iterator begin() const {
            return ConstIterator<T, binaryProxy>(*this, 0);
        }

        iterator end(){
            return Iterator<T, binaryProxy>(*this, size());
        }

        const_iterator end() const {
            return ConstIterator<T, binaryProxy>(*this, size());
        }

    private:
        const V1 &lhs;
        const V2 &rhs;
        const F  &func;
    };

    template<typename T, typename V, typename F>
    class binaryProxyRScalar{
    public:
        using value_type     = T;
        using iterator       = Iterator<T, binaryProxyRScalar>;
        using const_iterator = ConstIterator<T, binaryProxyRScalar>;

        binaryProxyRScalar(const V &left, Scalar<T> right, const F &f) 
            : lhs{left}, rhs{right}, func{f} {}

        T operator[](size_t i) const{
            return func(lhs[i], rhs[i]);
        }

        size_t size() const { return lhs.size(); }

        iterator begin(){
            return Iterator<T, binaryProxyRScalar>(*this, 0);
        }

        const_iterator begin() const{
            return ConstIterator<T, binaryProxyRScalar>(*this, 0);
        }

        iterator end(){
            return Iterator<T, binaryProxyRScalar>(*this, size());
        }

        const_iterator end() const{
            return ConstIterator<T, binaryProxyRScalar>(*this, size());
        }

    private:
        const V   &lhs;
        Scalar<T> rhs;
        const F   &func;
    };

    template<typename T, typename V, typename F>
    class binaryProxyLScalar{
    public:
        using value_type     = T;
        using iterator       = Iterator<T, binaryProxyLScalar>;
        using const_iterator = ConstIterator<T, binaryProxyLScalar>;

        binaryProxyLScalar(Scalar<T> left, const V &right, const F &f) 
            : lhs{left}, rhs{right}, func{f} {}

        T operator[](size_t i) const{
            return func(lhs[i], rhs[i]);
        }

        size_t size() const { return rhs.size(); }

        iterator begin(){
            return Iterator<T, binaryProxyLScalar>(*this, 0);
        }

        const_iterator begin() const{
            return ConstIterator<T, binaryProxyLScalar>(*this, 0);
        }

        iterator end(){
            return Iterator<T, binaryProxyLScalar>(*this, size());
        }

        const_iterator end() const{
            return ConstIterator<T, binaryProxyLScalar>(*this, size());
        }

    private:
        Scalar<T> lhs;    
        const V   &rhs;
        const F   &func;
    };

    /* unary operations */
    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, M, N, binaryProxyRScalar<T, V, std::plus<T>>> 
    operator+(const Matrix<T, M, N, V> &lhs, const T &elem){
        return binaryProxyRScalar<T, V, std::plus<T>>(lhs.data(), Scalar<T>{elem}, std::plus<T>());
    }

    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, M, N, binaryProxyLScalar<T, V, std::plus<T>>>
    operator+(const T &elem, const Matrix<T, M, N, V> &rhs){
        return binaryProxyLScalar<T, V, std::plus<T>>(Scalar<T>{elem}, rhs.data(), std::plus<T>());
    }

    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, M, N, binaryProxyRScalar<T, V, std::minus<T>>>
    operator-(const Matrix<T, M, N, V> &lhs, const T &elem){
        return binaryProxyRScalar<T, V, std::minus<T>>(lhs.data(), Scalar<T>{elem}, std::minus<T>());
    }

    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, M, N, binaryProxyLScalar<T, V, std::multiplies<T>>>
    operator*(const T &elem, const Matrix<T, M, N, V> &rhs){
        return binaryProxyLScalar<T, V, std::multiplies<T>>(Scalar<T>{elem}, rhs.data(), std::multiplies<T>());
    }

    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, M, N, binaryProxyRScalar<T, V, std::multiplies<T>>>
    operator*(const Matrix<T, M, N, V> &lhs, const T &elem){
        return binaryProxyRScalar<T, V, std::multiplies<T>>(lhs.data(), Scalar<T>{elem}, std::multiplies<T>());
    }

    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, M, N, binaryProxyRScalar<T, V, std::divides<T>>>
    operator/(const Matrix<T, M, N, V> &lhs, const T &elem){
        return binaryProxyRScalar<T, V, std::divides<T>>(lhs.data(), Scalar<T>{elem}, std::divides<T>());
    }    

    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, M, N, binaryProxyRScalar<T, V, std::modulus<T>>>
    operator%(const Matrix<T, M, N, V> &lhs, const T &elem){
        return binaryProxyRScalar<T, V, std::modulus<T>>(lhs.data(), Scalar<T>{elem}, std::modulus<T>());
    }       

    template<typename T, size_t M, size_t N, typename V1, typename V2>
    Matrix<bool, M, N, binaryProxy<bool, V1, V2, std::equal_to<T>>>
    operator==(const Matrix<T, M, N, V1> &lhs, const Matrix<T, M, N, V2> &rhs){
        return binaryProxy<bool, V1, V2, std::equal_to<T>>(lhs.data(), rhs.data(), std::equal_to<T>());
    }

    template<typename T, size_t M, size_t N, typename V1, typename V2>
    Matrix<bool, M, N, binaryProxy<bool, V1, V2, std::not_equal_to<T>>>
    operator!=(const Matrix<T, M, N, V1> &lhs, const Matrix<T, M, N, V2> &rhs){
        return binaryProxy<bool, V1, V2, std::not_equal_to<T>>(lhs.data(), rhs.data(), std::not_equal_to<T>());
    }

    template<typename T, size_t M, size_t N, typename V1, typename V2>
    Matrix<bool, M, N, binaryProxy<bool, V1, V2, std::greater<T>>>
    operator>(const Matrix<T, M, N, V1> &lhs, const Matrix<T, M, N, V2> &rhs){
        return binaryProxy<bool, V1, V2, std::greater<T>>(lhs.data(), rhs.data(), std::greater<T>());
    }

    template<typename T, size_t M, size_t N, typename V1, typename V2>
    Matrix<bool, M, N, binaryProxy<bool, V1, V2, std::greater_equal<T>>>
    operator>=(const Matrix<T, M, N, V1> &lhs, const Matrix<T, M, N, V2> &rhs){
        return binaryProxy<bool, V1, V2, std::greater_equal<T>>(lhs.data(), rhs.data(), std::greater_equal<T>());
    }

    template<typename T, size_t M, size_t N, typename V1, typename V2>
    Matrix<bool, M, N, binaryProxy<bool, V1, V2, std::less<T>>>
    operator<(const Matrix<T, M, N, V1> &lhs, const Matrix<T, M, N, V2> &rhs){
        return binaryProxy<bool, V1, V2, std::less<T>>(lhs.data(), rhs.data(), std::less<T>());
    }

    template<typename T, size_t M, size_t N, typename V1, typename V2>
    Matrix<bool, M, N, binaryProxy<bool, V1, V2, std::less_equal<T>>>
    operator<=(const Matrix<T, M, N, V1> &lhs, const Matrix<T, M, N, V2> &rhs){
        return binaryProxy<bool, V1, V2, std::less_equal<T>>(lhs.data(), rhs.data(), std::less_equal<T>());
    }

    template<typename T, size_t M, size_t N, size_t N1, typename V1, typename V2>
    class matrixMultiProxy{
    public:
        using value_type     = T;
        using iterator       = Iterator<T, matrixMultiProxy>;
        using const_iterator = ConstIterator<T, matrixMultiProxy>;

        matrixMultiProxy(const V1 &left, const V2 &right) 
            : lhs{left}, rhs{right} {}

        size_t size() const { return M*N1; }

        template<typename Q = V1>
        typename std::enable_if<MatrixImpl::IsParenType<Q>::value, T>::type
        operator()(size_t r, size_t c) const{
            T res = static_cast<T>(0);

            for(size_t i = 0; i < N; ++i)
                res += lhs(r, i)*rhs[i*N1+c];
            
            return res;
        }

        template<typename Q = V1>
        typename std::enable_if<!MatrixImpl::IsParenType<Q>::value, T>::type
        operator()(size_t r, size_t c) const{
            T res = static_cast<T>(0);

            for(size_t i = 0; i < N; ++i)
                res += lhs[r*N+i]*rhs[i*N1+c];
            
            return res;
        }        

        iterator begin(){
            return iterator(*this, 0);
        }

        const_iterator begin() const{
            return const_iterator(*this, 0);
        }

        iterator end(){
            return iterator(*this, size());
        }

        const_iterator end() const{
            return const_iterator(*this, size());
        }

    private:
        const V1 &lhs;
        const V2 &rhs;
    };

    template<typename T, size_t M1, size_t N, size_t N1, typename V1, typename V2>
    Matrix<T, M1, N1, matrixMultiProxy<T, M1, N, N1, V1, V2>>
    operator*(const Matrix<T, M1, N, V1> &lhs, const Matrix<T, N, N1, V2> &rhs){
        return matrixMultiProxy<T, M1, N, N1, V1, V2>(lhs.data(), rhs.data());
    }

    template<typename T, size_t M, typename V>
    class transProxy{
    public:
        using value_type     = T;
        using iterator       = Iterator<T, transProxy>;
        using const_iterator = ConstIterator<T, transProxy>;

        transProxy(const V &v) : vec{v} {}

        size_t size() const {
            return vec.size();
        }

        T operator()(size_t i, size_t j){
            return vec[j*M+i];
        }

        T operator()(size_t i, size_t j) const{
            return vec[j*M+i];
        }

        iterator begin(){
            return iterator(*this, 0);
        }

        const_iterator begin() const{
            return const_iterator(*this, 0);
        }

        iterator end(){
            return iterator(*this, size());
        }

        const_iterator end() const{
            return const_iterator(*this, size());
        }

    private:
        const V &vec;
    };

    // unary operations
    template<typename T, size_t M, size_t N, typename V>
    Matrix<T, N, M, transProxy<T, N, V>> 
    transpose(const Matrix<T, M, N, V> &m){
        return transProxy<T, N, V>(m.data());
    }
    
    template<typename T, size_t M, size_t N>
    class sliceMatrix{
    public:
        using iterator       = typename Matrix<T, M, N>::iterator;
        using const_iterator = typename Matrix<T, M, N>::const_iterator;

        sliceMatrix(Matrix<T, M, N> &m, slice row, slice col) 
            : mat{m}, r{row}, c{col} {}

        size_t rows() const { return r.size; }

        size_t cols() const { return c.size; }

        size_t size() const { return r.size*c.size; }

        T& operator()(size_t i, size_t j){
            MatrixImpl::index_bounds_check(*this, i, j);
            return mat.data()[r(i)+c.start+j*c.stride];            
        }

        const T& operator()(size_t i, size_t j) const{
            MatrixImpl::index_bounds_check(*this, i, j);
            return mat.data()[r(i)+c.start+j*c.stride];
        }

        sliceMatrix& operator=(std::initializer_list<T> il){
            assert(il.size() == size() && "overinput");

            std::vector<T> vec{il};
            for(size_t i = 0; i < rows(); ++i)
                for(size_t j = 0; j < cols(); ++j)
                    (*this)(i, j) = vec[i*cols()+j];

            return *this;
        }

        sliceMatrix& operator=(const sliceMatrix &rhs){
            for(size_t i = 0; i < rows(); ++i)
                for(size_t j = 0; j < cols(); ++j)
                    (*this)(i, j) = rhs(i, j);

            return *this;
        }

        iterator begin(){
            return mat.begin();
        }
        
        const_iterator begin() const{
            return mat.begin();
        }

        iterator end(){
            return mat.end();
        }

        const_iterator end() const{
            return mat.end();
        }

    private:
        Matrix<T, M, N> &mat;   
        slice r;
        slice c;
    };

    template<typename T, size_t M, size_t N>
    std::ostream& operator<<(std::ostream &os, const sliceMatrix<T, M, N> &sm){
        for(size_t i = 0; i < sm.rows(); ++i){
            for(size_t j = 0; j < sm.cols(); ++j){
                os.width(8);
                os.precision(3);
                os << std::fixed << sm(i, j) << " ";
            }
            os << '\n';
        }            
        os << '\n';
        return os;
    }

}   // Lee
