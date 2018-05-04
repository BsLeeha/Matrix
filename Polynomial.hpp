#ifndef _POLYNOMIAL
#define _POLYNOMIAL

#include <iostream>
#include <set>
#include <cmath>    // pow

namespace{
    struct Term{
        int coefficients;
        int exponent;

        Term(int c, int e) : coefficients{c}, exponent{e} {}
        // Term(std::initializer_list<int> il) 
        // : coefficients{*il.begin()}, exponent{*(il.begin()+1)}{};
    };

    struct TermComparator{
        bool operator()(const Term &lhs, const Term &rhs) const{
            return lhs.exponent < rhs.exponent;
        }
    };


    class Poly{
        friend std::ostream& operator<<(std::ostream &os, const Poly &p);
    public:
        Poly() : terms{}{};
        Poly(std::set<Term, TermComparator> &s) { terms = s; } 
        // Poly(std::initializer_list<std::initializer_list<int>> &il);
        ~Poly() {};

        Poly operator-() const;               // const: to operate on const object
        Poly operator+(const Poly &p);
        Poly operator-(const Poly &p);
        double operator()(double) const;
    private:
        std::set<Term, TermComparator> terms;
    };  

    std::ostream& operator<<(std::ostream &os, const Poly &p){
        int k = 1;
        os << "y = ";
        for(auto c : p.terms){
            if(k++!=1 && c.coefficients>0) os << "+";
            if(c.coefficients != 1 && c.coefficients != -1) os << c.coefficients;
            if(c.coefficients == -1) os << "-";
            if(c.exponent != 0) os << "x";
            if(c.exponent != 0 && c.exponent != 1) os << "^" << c.exponent;
        }
        os << "\n";
        return os;
    }

    // Poly::Poly(std::initializer_list<std::initializer_list<int>> &il){
    //     std::copy(il.begin(), il.end(), terms.begin());
    // }

    Poly Poly::operator-() const{
        Poly res;
        for(auto i : terms){
            res.terms.insert(Term(-i.coefficients, i.exponent));
        }
        return res;
    }

    Poly Poly::operator+(const Poly &p){
        auto my_it = terms.begin();
        auto p_it = p.terms.begin();
        Poly res;

        while(my_it != terms.end() && p_it != p.terms.end()){
            if(my_it->exponent < p_it->exponent){
                res.terms.insert(*my_it);
                ++my_it;
            }
            else if(my_it->exponent == p_it->exponent){
                res.terms.insert(Term(my_it->coefficients+p_it->coefficients, my_it->exponent));
                ++my_it;
                ++p_it;
            }
            else{
                res.terms.insert(*p_it);
                ++p_it;
            }
        }

        res.terms.insert(my_it, terms.end());
        res.terms.insert(p_it, p.terms.end());
        return res;
    }

    Poly Poly::operator-(const Poly &p){
        return (*this)+(-p);
    }

    double Poly::operator()(double t) const{
        double sum = 0;
        for(auto i : terms){
            sum += i.coefficients*pow(t, i.exponent);
        }
        return sum;
    }
}
#endif