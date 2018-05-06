#include "Matrix.hpp"
#include "Matrix_ref.hpp"

void Matrix_ref_test(){
    std::cout << "\nMatrix_ref test:\n";

    Lee::Matrix<double, 3> a{{{1, 2}, {3, 4}, {5, 6}}, {{7, 8}, {9, 10}, {11, 12}}};
    // Lee::Matrix_ref<a.description(), a.elem()> b;
    // std::cout << b;
    
}