#include <boost/rational.hpp>
#include <iostream>
#include <vector>
#include <boost/multiprecision/gmp.hpp>

using namespace std;

class Qn;
typedef boost::multiprecision::mpz_int integer;
typedef boost::rational<integer> fraction;
typedef vector<Qn> matrix;

namespace boost
{
    template <typename IntType>
    constexpr IntType round(rational<IntType> const &r)
    {
        return static_cast<IntType>((r.numerator() * 2 + r.denominator()) / (r.denominator() * 2));
    }
}

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

    void operator+=(const Qn &u)
    {
        assert(n == u.n && n != 0);
        for (int i = 0; i < n; i++)
            v[i] += u.v[i];
    }
    void operator-=(const Qn &u)
    {
        assert(n == u.n && n != 0);
        for (int i = 0; i < n; i++)
            v[i] -= u.v[i];
    }

    void operator=(const Qn &u)
    {
        assert(n == u.n && n != 0);
        for (int i = 0; i < n; i++)
            v[i] = u.v[i];
    }

    friend fraction dot(Qn &v, Qn &u);
    fraction normSqr() { return dot(*this, *this); }
};

Qn operator*(fraction f, Qn &u)
{
    Qn v(u);
    for (int i = 0; i < u.n; i++)
        v.v[i] *= f;
    return v;
}

fraction dot(Qn &v, Qn &u)
{
    fraction r(0);
    assert(v.n == u.n && v.n != 0);
    for (int i = 0; i < v.n; i++)
        r += v.v[i] * u.v[i];
    return r;
}

void init(matrix &x, int n, int m)
{
    assert(x.size() == 0);
    for (int i = 0; i < m; i++)
        x.push_back(Qn(n));
}

void swap(matrix &m, int k)
{
    fraction x;
    for (int i = 0; i < m[k].n; i++)
    {
        x = m[k].v[i];
        m[k].v[i] = m[k - 1].v[i];
        m[k - 1].v[i] = x;
    }
}

void MatrixIp(matrix &ret)
{
    ret.clear();
    int n, m;
    cout << "Enter nrows: ";
    cin >> n;
    cout << "Enter ncols: ";
    cin >> m;
    integer x;
    for (int i = 0; i < m; i++)
    {
        ret.push_back(Qn(n));
        for (int j = 0; j < n; j++)
        {
            cout << "Enter element [" << j << "," << i << "]: ";
            cin >> x;
            ret[i].v[j] = fraction(x);
        }
    }
}

void MatrixOp(matrix &ret, string s = "", int k = 0)
{
    int n = ret[0].n, m = ret.size();
    if (s != "")
        cout << s << endl;
    cout
        << "nrows: " << ret[0].n << endl;
    cout << "ncols: " << ret.size() << endl;
    integer x;
    cout << "[\n";
    for (int i = 0; i < n; i++)
    {
        cout << "[ ";
        for (int j = 0; j < m; j++)
        {
            if (k == 1)
                cout << ret[j].v[i] << " ";
            else
                cout << ret[j].v[i].numerator() << " ";
        }
        cout << "]\n";
    }
    cout << "]\n";
}

void GramSchmidt(matrix &B, matrix &mu, matrix &red) // returns B' : GS orthogonal basis
{
    int n = B.size();
    for (int i = 0; i < n; i++)
    {
        red[i] = (B[i]);
        mu[i].v[i] = fraction(integer(1));
        for (int j = 0; j < i; j++)
        {
            mu[i].v[j] = dot(B[i], red[j]) / red[j].normSqr();
            red[i] -= (mu[i].v[j] * red[j]);
        }
    }
}

void lll(matrix &basis, fraction delta)
{
    matrix mu, red;
    int n = basis.size(), k = 1;
    init(mu, n, n);
    init(red, n, n);
    GramSchmidt(basis, mu, red);

    while (k < n)
    {

        for (int j = k - 1; j >= 0; j--)
        {
            if (boost::round(mu[k].v[j]) != 0)
            {
                basis[k] -= fraction(boost::round(mu[k].v[j])) * basis[j];
                GramSchmidt(basis, mu, red);
            }
        }

        if (red[k].normSqr() >= (delta - (mu[k].v[k - 1] * mu[k].v[k - 1])) * (red[k - 1].normSqr()))
            k += 1;
        else
        {
            swap(basis, k);
            GramSchmidt(basis, mu, red);
            k = max(k - 1, 1);
        }
    }
    MatrixOp(basis, "Reduced Basis");
}

int main()
{
    matrix x;
    Qn t(4);
    t.v[0] = integer(105);
    t.v[1] = integer(821);
    t.v[2] = integer(404);
    t.v[3] = integer(328);
    x.push_back(t);
    t.v[0] = integer(881);
    t.v[1] = integer(667);
    t.v[2] = integer(644);
    t.v[3] = integer(927);
    x.push_back(t);
    t.v[0] = integer(181);
    t.v[1] = integer(483);
    t.v[2] = integer(87);
    t.v[3] = integer(500);
    x.push_back(t);
    t.v[0] = integer(893);
    t.v[1] = integer(834);
    t.v[2] = integer(732);
    t.v[3] = integer(441);
    x.push_back(t);
    lll(x, fraction(integer(3), integer(4)));
}

//[105, 821, 404, 328], [881, 667, 644, 927], [181, 483, 87, 500], [893, 834, 732, 441]