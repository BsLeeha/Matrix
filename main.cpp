#include <iostream>
#include <tuple>
#include <vector>
#include <array>
#include <set>
#include <complex>
#include <cmath>
#include <chrono>
#include <algorithm>
#include "Matrix_Test.cpp"
#include "Matrix_Desc_Test.cpp"
#include "Matrix_ref_test.cpp"

using namespace std;

int main(){
try{
    auto t0 = chrono::high_resolution_clock::now();

    Matrix_Desc_Test();
    Matrix_Test();
    Matrix_ref_test();

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