#include "Matrix.hpp"

void Matrix_Test(){
    std::cout << "\nMatrix Test:\n";
    
    // Lee::Matrix<int, 0> m1 = 5;
    // Lee::Matrix<int, 1> m2 {2, 3};
    Lee::Matrix<int, 2> m3 {{1, 2, 3}, {3, 4, 5}, {4, 5, 6}};
    Lee::Matrix<int, 2> m4 {{7, 8, 9}, {10, 11, 12}, {13, 14, 15}};
    // Lee::Matrix<int, 3> m5 {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};
    // Lee::Matrix<int, 3> m6 {{{9, 10}, {11, 12}}, {{13, 14}, {15, 16}}};    
    // cout << m5 << m6 << '\n';
    // cout << m5.row(0).col(0) << '\n';
    // cout << m5.col(0).row(0);

    Lee::Matrix<int, 2> m5(100, 100);
    srand(time(NULL));
    for(auto &i : *m5.elem()) i = rand()%100;

    Lee::Matrix<int, 2> m6(100, 100);
    srand(time(NULL));
    for(auto &i : *m6.elem()) i = rand()%100;

    std::cout << m5*m6; 
}