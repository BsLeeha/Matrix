#include "Matrix.hpp"
#include "Elimination.hpp"

void Matrix_Test(){
    std::cout << "\nMatrix Test:\n";
    
    // Lee::Matrix<int, 0> m1 = 5;
    // Lee::Matrix<int, 1> m2 {2, 3};
    Lee::Real_Matrix m3 {{1, 2, 3}, {3, 4, 5}, {4, 5, 6}};
    Lee::Real_Matrix m4 {{7, 8, 9}, {10, 11, 12}, {13, 14, 15}};
    std::cout << m3;
    m3.colPermute(1, 2);
    std::cout << m3;

    Lee::Matrix<int, 3> m5 {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};
    Lee::Matrix<int, 3> m6 {{{9, 10}, {11, 12}}, {{13, 14}, {15, 16}}};    
    // cout << m5 << m6 << '\n';
    // cout << m5.row(0).col(0) << '\n';        // attention
    // cout << m5.col(0).row(0);
    // std::cout << m5;
    // m5.row(0).col(0).rowPermute(0, 1);
    // std::cout << m5;

    // Lee::Matrix<float, 2> m5(500, 500);
    // srand(time(NULL));
    // for(auto &i : *m5.elem()) i = rand()%100;

    // Lee::Matrix<float, 2> m6(500, 500);
    // srand(time(NULL));
    // for(auto &i : *m6.elem()) i = rand()%100;

    // m5*m6; 
}