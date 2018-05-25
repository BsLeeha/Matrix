#include "Matrix.hpp"
#include "Matrix_Impl.hpp"

using namespace std;

void Matrix_Test(){
    std::cout << "\nMatrix Test:\n";
    
    // Lee::RealScalar m{2};
    // Lee::RealVectorr<3> m1 = {1, 2, 3};
    // m1 = {2, 3, 4};
    // std::cout << m1;
    // Lee::Matrix<double, 3, 3> m3 = {{1, 2, 3}, {3, 4, 5}, {4, 5, 6}};
    // Lee::Matrix<double, 3, 3> m4 {{7, 8, 9}, {10, 11, 12}, {13, 14, 15}};
    // std::cout << (m3 + (Lee::Matrix<double, 3, 3>{{1, 2, 3}, {4, 5, 6}, {4, 5, 6}}));
    const Lee::Matrix<int, 3, 2> m{{1, 2}, {3}};
    //m = {{1, 2}, {3, 4}};
    // cout << boolalpha << MatrixImpl::IsMatrixType<Lee::Matrix<int, 3, 2>>::value;
    MatrixImpl::test();

    // Lee::Matrix<int, 3> m5 {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};
    // Lee::Matrix<int, 3> m6 {{{9, 10}, {11, 12}}, {{13, 14}, {15, 16}}};    
    // cout << m5 << m6 << '\n';
    // cout << m5.row(0).col(0) << '\n';
    // cout << m5.col(0).row(0);

    // Lee::Matrix<int, 2> m5(100, 100);
    // srand(time(NULL));
    // for(auto &i : *m5.elem()) i = rand()%100;

    // Lee::Matrix<int, 2> m6(100, 100);
    // srand(time(NULL));
    // for(auto &i : *m6.elem()) i = rand()%100;
}