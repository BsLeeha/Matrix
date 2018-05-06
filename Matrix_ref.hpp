#ifndef _MATRIX_REF_H
#define _MATRIX_REF_H

#include "Matrix_Desc.hpp"

namespace Lee{
    template<typename T, size_t N>
    class Matrix_ref{     
    template<typename U>
    friend std::ostream& operator<<(std::ostream &os,  Matrix_ref<U, 1> m);        
    template<typename U, size_t M>
    friend std::ostream& operator<<(std::ostream &os,  Matrix_ref<U, M> m);    
    public:
        Matrix_ref(Lee::Matrix_desc<N> d, std::vector<T>* p)
            : desc{d}, elems{p}
        {}

        // basic info
        static constexpr size_t order() { return N; }
        size_t size() const { return desc.size; }
        size_t extention(size_t i) const { return desc.extents[i]; }
        Lee::Matrix_desc<N> description(){ return desc; }

        size_t rows() const { return extention(0); }
        size_t cols() const { return extention(1); }        

        // fortan style elem access
        template<typename...Dims>
        typename std::enable_if<std::is_convertible<Dims..., size_t>::value, T&>::type 
        operator()(Dims...dims){
            return (*elems)[desc(dims...)];
        }

        template<typename...Dims>
        typename std::enable_if<std::is_convertible<Dims..., size_t>::value, const T&>::type
        operator()(Dims...dims) const{
            return (*elems)[desc(dims...)];
        }

        // row/col access
        Matrix_ref<T, N-1>       operator[](size_t i) { return row(i); }
        Matrix_ref<T, N-1>       operator[](size_t i) const { return row(i); }

        Matrix_ref<T, N-1>       row(size_t i){
            static_assert(N-1>0, "1 dim matrix has no row to choose");
            assert(i < rows() && "row index over bound");

            Matrix_desc<N-1> slice;
            slice.start = desc.start + i*desc.strides[0];
            std::copy(desc.extents.begin()+1, desc.extents.end(), slice.extents.begin());
            slice.strides_init();

            Matrix_ref<T, N-1> res{slice, elems};
            return res;
        }

        Matrix_ref<const T, N-1> row(size_t i) const{
            static_assert(N-1>0, "1 dim matrix has no row to choose");            
            assert(i < rows() && "row index over bound");

            Matrix_desc<N-1> slice;
            slice.start = desc.start + i*desc.strides[0];
            std::copy(desc.extents.begin()+1, desc.extents.end(), slice.extents.begin());
            slice.strides_init();

            Matrix_ref<const T, N-1> res{slice, elems};
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

            Matrix_ref<T, N-1> res{coldesc, elems};
            return res;             
        }        

        Matrix_ref<T, N-1>       col(size_t i) const{
            static_assert(N-1>0, "1 dim matrix has no col to choose");            
            assert(i < cols() && "col index over bound");
       
            Matrix_desc<N-1> coldesc;
            coldesc.start = desc.start + i*desc.strides[1];
            std::copy(desc.extents.begin(), desc.extents.begin()+1, coldesc.extents.begin());
            if(N != 2) std::copy(desc.extents.begin()+2, desc.extents.end(), coldesc.extents.begin()+1);
            coldesc.size = desc.size/desc.extents[1];
            std::copy(desc.strides.begin(), desc.strides.begin()+1, coldesc.strides.begin());
            if(N != 2) std::copy(desc.strides.begin()+2, desc.strides.end(), coldesc.strides.begin()+1); 

            Matrix_ref<T, N-1> res{coldesc, elems};
            return res;             
        }         

        // row permute
        void rowPermute(size_t i, size_t j){
            for(size_t k = 0; k < desc.size/desc.extents[0]; ++k){
                T temp = (*elems)[desc.start + i*desc.strides[0] + k];
                
                (*elems)[desc.start + i*desc.strides[0] + k] = 
                (*elems)[desc.start + j*desc.strides[0] + k];
                
                (*elems)[desc.start + j*desc.strides[0] + k] = temp;
            }
        }

        // col permute
        void colPermute(size_t i, size_t j){
            for(size_t r = 0; r < rows(); ++r){
                for(size_t k = 0; k < desc.size/desc.extents[0]/desc.extents[1]; ++k){
                    T temp = (*elems)[desc.start + r*desc.strides[0] + i*desc.strides[1] + k];

                    (*elems)[desc.start + r*desc.strides[0] + i*desc.strides[1] + k] = 
                    (*elems)[desc.start + r*desc.strides[0] + j*desc.strides[1] + k];

                    (*elems)[desc.start + r*desc.strides[0] + j*desc.strides[1] + k] = temp;
                }
            }
        }

        // value change
        void operator=(const Matrix_ref<T, 1> &rhs){
            assert(desc.extents == rhs.desc.extents && "extents not match");

            for(size_t i = 0; i < size(); ++i){
                (*this)(i) = rhs(i);
            }
        }   

        template<size_t U>
        void operator=(const Matrix_ref<T, U> &rhs){
            assert(desc.extents == rhs.desc.extents && "extents not match");

            for(size_t i = 0; i < rows(); ++i){
                (*this)[i] = rhs[i];
            }            
        }

    private:
        Lee::Matrix_desc<N> desc;
        std::vector<T>* elems;
    };

    template<typename T>
    std::ostream& operator<<(std::ostream &os,  Matrix_ref<T, 1> m){
        os << '{';
        for(size_t i = 0; i < m.size(); ++i)
        {
            os << m(i);
            if(i != m.size()-1) os << ", ";
        }
        os << '}';
        return os;
    }  

    template<typename T, size_t N>
    std::ostream& operator<<(std::ostream &os,  Matrix_ref<T, N> m){
        os << '{';
        for(size_t i = 0; i < m.rows(); ++i){
            os << m[i];
            if(i+1 != m.rows()) os << ", ";
        }
        os << "}";
        return os;
    }

    // template<typename T>
    //     void operator=(Matrix_ref<T, 1> &lhs, const Matrix_ref<T, 1> &rhs){
    //         assert(lhs.desc.extents == rhs.desc.extents && "extents not match");

    //         for(size_t i = 0; i < size(); ++i){
    //             lhs(i) = rhs(i);
    //         }
    //     }   

    // template<typename T, size_t N>
    //     void operator=(Matrix_ref<T, N> &lhs, const Matrix_ref<T, N> &rhs){
    //         assert(lhs.desc.extents == rhs.desc.extents && "extents not match");
            
    //         for(size_t i = 0; i < rows(); ++i){
    //             lhs[i] = rhs[i];
    //         }            
    //     }

}   // Lee

#endif  // _MATRIX_REF_H