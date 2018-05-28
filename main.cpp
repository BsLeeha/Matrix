#include <iostream>
#include <tuple>
#include <vector>
#include <array>
#include <set>
#include <complex>
#include <cmath>
#include <chrono>
#include <algorithm>
#include "Matrix.hpp"
#include "Matrix_Impl.hpp"

using namespace std;

int main(){
try{
    auto t0 = chrono::high_resolution_clock::now();

    Lee::Matrix<int, 3, 2> m(6, -1);
    Lee::Matrix<int, 3, 2> m1{1, 2, 3, 4};
    Lee::Matrix<int, 2, 3> m5{1, 2, 3, 4};
    Lee::Matrix<int, 3, 2> m2{0, 0, 4, 3, 2, 1};
    Lee::Matrix<int, 3, 2> m3 = 2*(1+2+m-1+2+3)*2/16%2+10;
    Lee::Matrix<int, 3, 2> m4 = m1(Lee::slice(0, 2, 2), Lee::slice(0, 1, 1));
    std::cout << m << Lee::abs(m) << m3;

    auto t1 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
    cout << "\nduration time: " << duration << " milliseconds\n";
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
