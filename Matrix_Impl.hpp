#pragma once

#include <type_traits>
#include <algorithm>
#include <functional>

namespace Lee{
    template<typename T, size_t M, size_t N, typename V>
    class Matrix;    

    template<typename T, size_t M, typename V>
    class transProxy;    

    template<typename T, size_t M, size_t N, size_t N1, typename V1, typename V2>
    class matrixMultiProxy;    

}

namespace MatrixImpl{

    template<typename T, size_t N>
    struct NestedInitializerList{
        using type = std::initializer_list<typename NestedInitializerList<T, N-1>::type>;
    };

    template<typename T>
    struct NestedInitializerList<T, 1>{
        using type = std::initializer_list<T>;
    };    

    template<typename T>
    struct NestedInitializerList<T, 0>;        // no such initializer_list

    template<typename T, size_t N>
    using NestedInitializerListN = typename NestedInitializerList<T, N>::type;

    template<typename T>
    struct IsMatrixType{
        static const bool value = false;
    };

    template<typename T, size_t M, size_t N, typename V>
    struct IsMatrixType<Lee::Matrix<T, M, N, V>>{
        static const bool value = true;
    }; 

    template<typename T>
    struct IsParenType{
        static const bool value = false;
    };

    template<typename T, size_t M, typename V>
    struct IsParenType<Lee::transProxy<T, M, V>>{
        static const bool value = true;
    };

    template<typename T, size_t M, size_t N, size_t N1, typename V1, typename V2>
    struct IsParenType<Lee::matrixMultiProxy<T, M, N, N1, V1, V2>>{
        static const bool value = true;
    };



    template<typename M>
    void index_bounds_check(const M &m, size_t r, size_t c){
        if(IsMatrixType<M>::value) assert(r<m.rows() && c<m.cols() && "index out of range");
    }    

    template<typename M>
    void matrix_valid(const M &m){
        assert(m.data().size() == m.rows()*m.cols() && "matrix is not valid");
    }

    template<typename M, typename...Dims>
    void matrix_valid(const M &m, const Dims &...dims){
        assert(m.data().size() == m.rows()*m.cols() && "matrix is not valid");
        matrix_valid(dims...);
    }


    constexpr bool Request_index() { return true; }

    template<typename T, typename...Args>
    constexpr bool Request_index(T first, Args...args){
        return std::is_convertible<T, size_t>::value && Request_index(args...);
    }

    template<typename T, typename...Dims>
    bool check_bounds(T slice, Dims...dims){
        size_t N = sizeof...(Dims);
        size_t arr[N]{size_t(dims)...};
        for(size_t i = 0; i < N; ++i)
            if(arr[i] >= slice.extents[i]) return false;
        return true;
    }

}   //Matrix_Impl

