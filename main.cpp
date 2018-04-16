#include <iostream>
#include <tuple>
#include <set>
#include <cmath>
#include "Matrix.hpp"
#include "Elimination.hpp"
#include "Factorization.hpp"
#include "EigenvalueEstimate.hpp"
#include "Polynomial.hpp"
#include "SystemSolving.hpp"

using namespace std;

int main(){
try{
    Lee::Matrix<double, 3, 3> e{3, 1, -1, 2, 4, 1, -1, 2, 5};
    Lee::Matrix<double, 3, 4> e1{1, 2, 2, 2, 2, 4, 6, 8, 3, 6, 8, 10};
    Lee::Matrix<double, 3, 1> b{4, 1, 1};
    Lee::Matrix<double, 3, 2> A = {1, -4, 2, 3, 2, 2};
    Lee::Matrix<double, 2, 1> x = {0, 1};
    vector<double> pivots = pivot(A);
    cout << A << "\n" << get<0>(QRx(A)) << "\n" << get<1>(QRx(A)) << "\n" << get<0>(QRx(A))*get<1>(QRx(A));
    // cout << DJI(A, b, Lee::Matrix<double, 2, 1>{1, 1}, 0.0005);
    // set<Term, TermComparator> s1 = {{1, 3}, {1, 1}, {-1, 0}};
    // Poly p1{s1};
    // cout << bisect(p1, 0, 1, 0.0005);
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