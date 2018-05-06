#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <array>
#include <cstdlib>
#include <tuple>
#include <algorithm>
#include <type_traits>
#include <cassert>
#include "Matrix_Impl.hpp"
#include "Matrix_Desc.hpp"
#include "Matrix_ref.hpp"

namespace Lee{

    template<typename T, size_t N>
    class Matrix_base{
    public:
    private:
    
    };

    template<typename T, size_t N>
    using Matrix_initializer = MatrixImpl::NestedInitializerListN<T, N>;

    template<typename T, size_t N>      // N dimension
    class Matrix{
    template<typename U>
    friend std::ostream& operator<<(std::ostream &os,  Matrix<U, 1> m);
    template<typename U, size_t M>
    friend std::ostream& operator<<(std::ostream &os, Matrix<U, M> m);

    public:
        // Conventional aliases
        using value_type = T;
        using iterator = typename std::vector<T>::iterator;
        using const_iterator = typename std::vector<T>::const_iterator;
        using reverse_iterator = typename std::vector<T>::reverse_iterator;
        using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

        // Initialization
        Matrix() = default;                                         // default
        Matrix(const Matrix &) = default;                           // copy
        Matrix& operator=(const Matrix&) = default;     
        Matrix(Matrix &&) = default;                                // move
        Matrix& operator=(Matrix&&) = default;
        ~Matrix() = default;

        template<typename...Ts>                                     // (), initialize dims from extents
        Matrix(Ts...ts)
            : desc{0, ts...}, elems(MatrixImpl::product(ts...)){}

        Matrix(Matrix_initializer<T, N> list){                      // {}, initialize elems from nested initializer_list
            desc.start = 0;
            MatrixImpl::derive_extents<N>(list, desc.extents.data());
            desc.strides_init();
            elems.reserve(desc.size);
            MatrixImpl::derive_elems<T, N>(list, elems);

            assert(desc.size == elems.size() && "extents size and elems size not match");
        }
        Matrix& operator=(Matrix_initializer<T, N>);

        // template<typename U, typename std::enable_if<(N>1), void>::type* = nullptr>
        // Matrix(std::initializer_list<U>) = delete;                  // disable {} initialization with extents

        // template<typename U>
        // Matrix& operator=(std::initializer_list<U>) = delete;

        // basic info
        static constexpr size_t order() { return N; }
        const Matrix_desc<N>& description() const { return desc; }
        size_t extention(size_t i) const { return desc.extents[i]; }

        size_t size() const { return desc.size; }
        size_t rows() const { return extention(0); }
        size_t cols() const { return extention(1); }

        std::vector<T>* elem() { return &elems; }
        T*       data() { return elems.data(); }
        const T* data() const { return elems.data(); }    

        iterator begin() { return elems.begin(); }
        const_iterator begin() const { return elems.begin(); }
        const_iterator cbegin() const {return elems.cbegin(); }
        
        iterator end() { return elems.end(); }
        const_iterator end() const { return elems.end(); }        
        const_iterator cend() const { return elems.cend(); }            


        // fortan style elem access
        template<typename...Dims>
        T& operator()(Dims...dims){
            assert(MatrixImpl::Request_index(dims...) && "index type wrong");

            return elems[desc(dims...)];
        }

        template<typename...Dims>
        const T& operator()(Dims...dims) const{
            assert(MatrixImpl::Request_index(dims...) && "index type wrong");
            
            return elems[desc(dims...)];
        }

        // row/col access
        Matrix_ref<T, N-1>       operator[](size_t i) { return row(i); }
        Matrix_ref<const T, N-1> operator[](size_t i) const { return row(i); }

        Matrix_ref<T, N-1>       row(size_t i){
            static_assert(N-1>0, "1 dim matrix has no row to choose");
            assert(i < rows() && "row index over bound");

            Matrix_desc<N-1> slice;
            slice.start = desc.start + i*desc.strides[0];
            std::copy(desc.extents.begin()+1, desc.extents.end(), slice.extents.begin());
            slice.strides_init();

            Matrix_ref<T, N-1> res{slice, &elems};
            return res;
        }

         Matrix_ref<const T, N-1> row(size_t i) const{
            static_assert(N-1>0, "1 dim matrix has no row to choose");            
            assert(i < rows() && "row index over bound");

            Matrix_desc<N-1> slice;
            slice.start = desc.start + i*desc.strides[0];
            std::copy(desc.extents.begin()+1, desc.extents.end(), slice.extents.begin());
            slice.strides_init();

            Matrix_ref<T, N-1> res{slice, &elems};
            return res;        
        }

        Matrix_ref<T, N-1>       col(size_t i){
            static_assert(N-1>0, "1 dim matrix has no col to choose");            
            assert(i < cols() && "col index over bound");
       
            Matrix_desc<N-1> coldesc;
            coldesc.start = desc.start + i*desc.strides[1];
            std::copy(desc.extents.begin(), desc.extents.begin()+1, coldesc.extents.begin());
            if(N != 2) std::copy(desc.extents.begin()+2, desc.extents.end(), coldesc.extents.begin()+1);
            coldesc.size = desc.size/desc.extents[1];
            std::copy(desc.strides.begin(), desc.strides.begin()+1, coldesc.strides.begin());
            if(N != 2) std::copy(desc.strides.begin()+2, desc.strides.end(), coldesc.strides.begin()+1);            

            Matrix_ref<T, N-1> res{coldesc, &elems};
            return res;             
        }

        Matrix_ref<const T, N-1> col(size_t i) const{
            static_assert(N-1>0, "1 dim matrix has no col to choose");            
            assert(i < cols() && "col index over bound");
       
            Matrix_desc<N-1> coldesc;
            coldesc.start = desc.start + i*desc.strides[1];
            std::copy(desc.extents.begin(), desc.extents.begin()+1, coldesc.extents.begin());
            if(N != 2) std::copy(desc.extents.begin()+2, desc.extents.end(), coldesc.extents.begin()+1);
            coldesc.size = desc.size/desc.extents[1];
            std::copy(desc.strides.begin(), desc.strides.begin()+1, coldesc.strides.begin());
            if(N != 2) std::copy(desc.strides.begin()+2, desc.strides.end(), coldesc.strides.begin()+1);            

            Matrix_ref<T, N-1> res{coldesc, &elems};
            return res;             
        }

        // matrix equal judge
        bool operator==(const Matrix<T, N> &m) const{
            MatrixImpl::assert_size_equal(*this, m);

            return *data() == *m.data();
        }

        // matrix add
        Matrix<T, N>& operator+=(const Matrix<T, N> &m){
            MatrixImpl::assert_size_equal(*this, m);

            std::transform(begin(), end(), m.begin(), elems.begin(), std::plus<T>());
            return *this;
        }

        // matrix substract
        Matrix<T, N>& operator-=(const Matrix<T, N> &m){
            MatrixImpl::assert_size_equal(*this, m);

            std::transform(begin(), end(), m.begin(), elems.begin(), std::minus<T>());
            return *this;
        }

        // scalar multiplication
        Matrix<T, N>& operator*=(T m){
            std::for_each(begin(), end(), [m](T &item){ item *= m;});
            return *this;
        }

        // matrix multiplication
        Matrix<T, 2> operator*=(const Matrix<T, 2> &m){
            assert(cols() == m.rows() && "cannot multiply");

            Matrix<T, 2> res(rows(), m.cols());
            for(size_t i = 0; i < rows(); ++i)
                for(size_t j = 0; j < m.cols(); ++j)
                    for(size_t k = 0; k < cols(); ++k){
                        res(i, j) += (*this)(i, k) * m(k, j);
                    }
            return res;
        }

        // row permute
        void rowPermute(size_t i, size_t j){
            for(size_t k = 0; k < desc.size/desc.extents[0]; ++k){
                T temp = elems[desc.start + i*desc.strides[0] + k];
                
                elems[desc.start + i*desc.strides[0] + k] = 
                elems[desc.start + j*desc.strides[0] + k];
                
                elems[desc.start + j*desc.strides[0] + k] = temp;
            }
        }

        // col permute
        void colPermute(size_t i, size_t j){
            for(size_t r = 0; r < rows(); ++r){
                for(size_t k = 0; k < desc.size/desc.extents[0]/desc.extents[1]; ++k){
                    T temp = elems[desc.start + r*desc.strides[0] + i*desc.strides[1] + k];

                    elems[desc.start + r*desc.strides[0] + i*desc.strides[1] + k] = 
                    elems[desc.start + r*desc.strides[0] + j*desc.strides[1] + k];

                    elems[desc.start + r*desc.strides[0] + j*desc.strides[1] + k] = temp;
                }
            }
        }

        void to_eye() { std::for_each(elems.begin(), elems.end(), [](T &item){ item = 1; }); }
        void to_zero() { std::for_each(elems.begin(), elems.end(), [](T &item){ item = 0; }); }

    private:
        Lee::Matrix_desc<N> desc;                // slice descripting extents of each dimension    
        std::vector<T> elems;                     // elements
    };

    // overload ouput
    template<typename T>
    std::ostream& operator<<(std::ostream &os,  Matrix<T, 1> m){
        os << m.desc;        
        os << '{';
        for(size_t i = 0; i < m.size(); ++i){
            os << m(i);
            if(i+1 != m.size()) os << ", ";
        }
        os << "}\n";
        return os;
    }  

    template<typename T, size_t N>
    std::ostream& operator<<(std::ostream &os,  Matrix<T, N> m){
        os << m.desc;
        for(size_t i = 0; i < m.rows(); ++i){
            os << m[i] << '\n';
        }
        return os;
    }       

    /* Scalar Class Specialization */
    template<typename T>
    class Matrix<T, 0>{
    public:
        using value_type = T;

        // initialize
        Matrix() = default;
        Matrix(const Matrix<T, 0>&) = default;
        Matrix& operator=(const Matrix<T, 0> &) = default;
        Matrix(const T &x) : elem{x} {}
        Matrix& operator=(const T &x) { elem = x; return *this; }

        // elem access
        T& operator()() { return elem; }
        const T& operator()() const { return elem; }

    private:
        T elem;  
    };

    template<typename T>
    std::ostream& operator<<(std::ostream& os, const Matrix<T, 0> &m){
        os << m();
        return os;
    }

    /* Random Class for test */
    template<typename T, size_t N>
    class RandomMatrix : public Matrix<T, N>{
    public:
        RandomMatrix(size_t m, size_t n)
        : Matrix<T, N>(m, n){
            srand(time(NULL));
            for(auto &i : *Matrix<T, N>::elem()) i = rand()%100;
        }     
    };

    /* Convenient Alias */
    using Real_Scalar = Matrix<double, 0>;

    using Real_Vector = Matrix<double, 1>;

    using Real_Matrix = Matrix<double, 2>;


    // Matrix add
    template<typename T, size_t N>
    Matrix<T, N> operator+(Matrix<T, N> lhs, const Matrix<T, N> &rhs){
        return lhs += rhs;
    }

    // Matrix subtract
    template<typename T, size_t N>
    Matrix<T, N> operator-(Matrix<T, N> lhs, const Matrix<T, N> &rhs){
        return lhs -= rhs;
    }

    // Matrix multiplication
    template<typename T, size_t N>
    Matrix<T, N> operator*(Matrix<T, N> lhs, const Matrix<T, N> &rhs){
        return lhs *= rhs;
    }    

}   // Lee

#endif  // MATRIX_H