#ifndef _MATRIX_SLICE_H
#define _MATRIX_SLICE_H

#include "Matrix_Impl.hpp"
#include <numeric>
#include <type_traits>

/*
*  Matrix_slice: define the shape of a matrix(extents) and the mapping in memory(strides)
*/

namespace Lee{ 

    template<size_t N>
    struct Matrix_slice{

        // default initialization
        Matrix_slice()
            : start{0}, size{0}
            {
                extents.fill(0);
                strides.fill(0);
            }

        // initialize with extents provided
        Matrix_slice(size_t s, std::initializer_list<size_t> exts)
            : start{s}
        {
            assert(exts.size() == N && "extents size not meet");    // initializer_list.size() is constexpr in C++14 

            std::copy(exts.begin(), exts.end(), extents.begin());
            strides_init();               
        }

        template<typename...Exts>               
            Matrix_slice(size_t s, Exts...exts)
                : start{s}
            {
                static_assert(sizeof...(Exts) == N, "extents size not meet");           

                extents = {size_t(exts)...};
                strides_init();                
            }
        

        // extents and strides provided
        Matrix_slice(size_t s, std::initializer_list<size_t> exts, std::initializer_list<size_t> strs)
            : start{s}
        {
            assert(exts.size() == N && "extents size not meet");

            extents = exts;
            strides = strs;
            size = extents[0] * strides[0];
        } 

        // subscript: dims to index in memory

        template<typename...Dims>
        size_t operator()(Dims...dims) {
            static_assert(sizeof...(Dims) == N, "number of index not match");
            assert(MatrixImpl::Request_index(dims...) && "index type wrong");
            assert(MatrixImpl::check_bounds(*this, dims...) && "index over bounds");   

            size_t arr[N] {size_t(dims)...};
            return start + std::inner_product(arr, arr+N, strides.begin(), size_t(0));            
        }

        template<typename...Dims>
        size_t operator()(Dims...dims) const{
            static_assert(sizeof...(Dims) == N, "number of index not match");
            assert(MatrixImpl::Request_index(dims...) && "index type wrong");    
            assert(MatrixImpl::check_bounds(*this, dims...) && "index over bounds");

            size_t arr[N] {size_t(dims)...};
            return start + std::inner_product(arr, arr+N, strides.begin(), size_t(0));            
        }

        // size and strides deduction by extents
        void strides_init(){
            strides[N-1] = 1;
            for(size_t i = N-1; i > 0; --i){
                strides[i-1] = strides[i] * extents[i];
            }
            size = extents[0]*strides[0];            
        }

        size_t start;                           // starting offset
        size_t size;                            // total size
        std::array<size_t, N> extents;          // # elements in each dimension
        std::array<size_t, N> strides;          // offsets between elems in each dimension
    };

    template<size_t N>
    std::ostream& operator<<(std::ostream &os, const Matrix_slice<N> &m){
        os << "size: " << std::dec << m.size;
        os << " start: " << m.start; 
        os << " extents: ";
        for(auto i : m.extents) os << i << ' ';
        os << " strides: ";
        for(auto i : m.strides) os << i << ' ';
        os << '\n';
        return os;
    }

    // template<> template<>
    // size_t Matrix_slice<1>::operator()<>(size_t i)const{
    //     std::cout << "yes";     // specialization choosing fail
    //     return i;
    // }

    // template<> template<>
    // size_t Matrix_slice<2>::operator()<>(size_t i, size_t j)const{
    //     std::cout << "successful";  // specialization choosing fail
    //     return i*strides[0]+j;
    // }

    
}   // Lee

#endif  // _MATRIX_SLICE_H