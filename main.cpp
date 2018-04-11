#include <iostream>
#include <tuple>
#include "Matrix.hpp"
#include "Elimination.hpp"
#include "Factorization.hpp"

using namespace std;

int main(){
try{
    Lee::Matrix<double, 3, 3> e{1, 3, 0, 2, 1, 0, 1, 2, 0};
    Lee::Matrix<double, 3, 4> e1{1, 2, 2, 2, 2, 4, 6, 8, 3, 6, 8, 10};
    Lee::Matrix<double, 3, 3> e2{1, 2, 3, 4, 5, 6, 7, 8, 9};
    Lee::Matrix<double, 3, 1> b{0, 0, 0};
    Lee::Matrix<double, 3, 3> A = e;
    // Lee::Matrix<double, 3, 3> R = rref(A);
    std::tuple<Lee::Matrix<double, 3, 3>, Lee::Matrix<double, 3, 3>> res = PLU(A);

    cout << "A:\n" << A << "\n";    
    cout << "Lower: \n" << std::get<0>(res) << "\n";    
    cout << "Upper: \n" << std::get<1>(res) << "\n";
    cout << "Make up: \n" << std::get<0>(res)*std::get<1>(res) << "\n";
    return 0;
}
catch(out_of_range){
    cerr << "range error!\n";
}
catch(length_error){
    cerr << "length error!\n";
}
catch(bad_alloc){
    cerr << "bad alloc!\n";
}
catch(...){
    cerr<<"unknown exception error\n";
}
    return 0;
}