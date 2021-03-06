#include <boost/rational.hpp>
#include <iostream>
#include <vector>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

class Qn;
typedef boost::multiprecision::mpz_int integer;
typedef boost::rational<integer> fraction;
typedef boost::multiprecision::mpf_float floating;
typedef std::vector<integer> Zn;
typedef std::vector<Qn> matrix;

class Qn
{
public:
    fraction *v = NULL;
    int n = 0;

public:
    Qn(int t)
    {
        n = t;
        v = new fraction[t]();
    }

    Qn(const Qn &u)
    {
        n = u.n;
        v = new fraction[n];
        for (int i = 0; i < n; i++)
            v[i] = u.v[i];
    }

    ~Qn() { delete[](v); }

    void operator-=(const Qn &u);
    void operator=(const Qn &u);
    friend fraction dot(Qn &v, Qn &u);
    fraction normSqr() { return dot(*this, *this); }
};

Qn operator*(fraction f, Qn &u);
fraction dot(Qn &v, Qn &u);
void QnIp(Qn &v, int n);

void init(matrix &x, int n, int m);
void MatrixIp(matrix &ret);
void MatrixOp(matrix &ret, std::string s = "", int k = 0);

void GramSchmidt(matrix &B, matrix &mu, matrix &red);
void lll(matrix &basis, fraction delta);