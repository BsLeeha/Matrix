#include <vector>

namespace Lee{
    template<typename T, size_t N>
    class Matrix{
    public:
        static constexpr size_t order = N;                      // dimension
        using value_type = T;
        using iterator = typename std::vector<T>::iterator;
        using const_iterator = typename std::vector<T>::const_iterator;

        Matrix() = default;
        Matrix(Matrix &&) = default;                    // move
        Matrix& operator=(Matrix &&) = default;
        Matrix(Matrix const&) = default;                // copy
        Matrix& operator=(Matrix const&) = default;
        ~Matrix() = default;

        template<typename U>
            Matrix<const Matrix_ref<U, N>&>;            // construct from Matrix_ref
        template<typename U>
            Matrix& operator=(const Matrix_ref<U, N>&); // assign from Matrix_ref
        
        template<typename...Exts>                       // specify the extents               
            explicit Matrix(Exts...exts);
        
        Matrix(Matrix_initializer<T, N>);               // initialize from list
        Matrix& operator=(Matrix_initializer<T, N>);    // assign from list

        template<typename U>
            Matrix<initializer_list<U>> = delete;        // don't use {} except for elements
        template<typename U>
            Matrix& operator=(initializer_list<U>) = delete;

        static constexpr size_t order() { return N; }           // number of dimensions
        size_t extent(size_t n) const { return desc.extents[n]; }   // #elements in the nth dimension
        size_t size() const { return elems.size(); }                // total number of elements
        const Matrix_slice<N>& descriptor() const { return desc; }  // the slice defining subscripting

        T* data(){ return elems.data(); }                       // "flat" element access
        const T* data() { return elems.data(); }

    private:
        Matrix_slice<N> desc;                   // slice descriptor defining extents in the N dimensions
        vector<T> elems;                       
    };

    template<typename T, size_t N>
        template<typename...Exts>
        Matrix<T, N>::Matrix(Exts...exts)
            : desc{exts...},        // copy extents        
            elems{desc.size}        // allocate desc.size elements and default initialize them
        {}

    template<typename T, size_t N>
        Matrix<T, N>::Matrix(Matrix_initializer<T, N> init){    
            Matrix_impl::derive_extents(init, desc.extents);    // deduce extents from initializer list
            elems.reserve(desc.size);                           // make room for slices
            Matrix_impl::insert_flat(init, elems);              // initialize from initializer list
            assert(elems.size() == desc.size);
        }

    template<typename T, size_t N>
        Matrix<T, N>::Matrix(const Matrix_ref<U, N>& x)
            : desc{x.desc}, elems{x.begin(), x.end()}
        {
            static_assert(Convertible<U, T>(), "Matrix constructor: incompatible element types");
        }
    
    template<typename T, size_t N>
        template<typename U>
        Matrix<T, N>& Matrix<T, N>::operator=(const Matrix_ref<U, N>& x){
            static_assert(Convertible<U, T>(), "Matrix =: incompatible element types");

            desc = x.desc;
            elems.assign(x.begin(), x.end());
            return *this;
        }
    
    template<size_t N>
    struct Matrix_slice{
        Matrix_slice() = default;

        Matrix_slice(size_t s, initializer_list<size_t> exts);
        Matrix_slice(size_t s, initializer_list<size_t> exts, initializer_list<size_t> strs);

        template<typename...Dims>
            Matrix_slice(Dimes...dims);

        template<typename...Dims, 
                 typename = Enable_if<All(Convertible<Dims, size_t>()...)>>
            size_t operator()(Dims...dims) const;
        
        size_t size;
        size_t start;
        array<size_t, N> extents;   // number of elemts
        array<size_t, N> strides;   
    }
}