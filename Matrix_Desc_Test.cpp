#include "Matrix_Desc.hpp"

void Matrix_Desc_Test(){
    std::cout << "\nMatrix Desc Test:\n";
    Lee::Matrix_desc<1> m{0, 3};
    Lee::Matrix_desc<3> m1{0, 3, 4, 5};
    Lee::Matrix_desc<2> m2(0, 3, 4);    
    std::cout << m << '\n' << m1 << '\n' << m2 << '\n';    
    std::cout << m(2);
    std::cout << m2(2, 1);
}