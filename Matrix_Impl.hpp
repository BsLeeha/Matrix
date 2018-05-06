#ifndef _MATRiX_IMPL_H
#define _MATRiX_IMPL_H

#include <type_traits>
#include <algorithm>
#include <functional>

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

    constexpr bool Request_index() { return true; }

    template<typename T, typename...Args>
    constexpr bool Request_index(T first, Args...args){
        return std::is_convertible<T, size_t>::value && Request_index(args...);
    }

    size_t product(size_t a){
        return a;
    }

    template<typename...Ts>
    size_t product(size_t a, Ts...ts){
        assert(Request_index(ts...) && "index type wrong"); 

        return a * product(ts...);
    }

    template<size_t N, typename List>
    void derive_extents(const List &list, size_t *p, typename std::enable_if<(N==1), void>::type* = 0){
        *p = list.size();
    }

    template<size_t N, typename List>
    void derive_extents(const List &list, size_t *p, typename std::enable_if<(N>1), void>::type* =0){
        *p = list.size();
        derive_extents<N-1>(*list.begin(), ++p);
    }    

    // template<typename T, size_t N, typename List>
    // void derive_elems(const List &list, T *p, typename std::enable_if<(N==1), void>::type* = 0){
    //     for(auto i : list){
    //         *p++ = i;
    //     }
    // }

    // template<typename T, size_t N, typename List>
    // void derive_elems(const List &list, T *p, typename std::enable_if<(N>1), void>::type* = 0){
    //     for(auto i = list.begin(); i != list.end(); ++i){
    //         derive_elems<T, N-1>(*i, p);
    //         if(i != list.end()-1) p += i->size();
    //     }
    // }

    template<typename T, size_t N, typename List, typename Vec>
    void derive_elems(const List &list, Vec &v, typename std::enable_if<(N==1), void>::type* = 0){
        for(auto i : list){
            v.push_back(i);
        }
    }

    template<typename T, size_t N, typename List, typename Vec>
    void derive_elems(const List &list, Vec &v, typename std::enable_if<(N>1), void>::type* = 0){
        for(auto i = list.begin(); i != list.end(); ++i){
            derive_elems<T, N-1>(*i, v);
        }
    }    

    template<typename T, typename...Dims>
    bool check_bounds(T slice, Dims...dims){
        size_t N = sizeof...(Dims);
        size_t arr[N]{size_t(dims)...};
        for(size_t i = 0; i < N; ++i)
            if(arr[i] >= slice.extents[i]) return false;
        return true;
    }

    template<typename M>
    void assert_size_equal(const M &lhs, const M &rhs) {
            assert((lhs.description() == rhs.description()) && "size not match");        
    }

}   //Matrix_Impl

#endif // _MATRiX_IMPL_H