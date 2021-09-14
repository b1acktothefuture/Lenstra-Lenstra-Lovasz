#include "lll.h"
namespace boost
{
    template <typename IntType>
    constexpr IntType round(rational<IntType> const &r)
    {
        return static_cast<IntType>((r.numerator() * 2 + r.denominator()) / (r.denominator() * 2));
    }
}

/****************************** Qn **********************************/

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

void Qn ::operator-=(const Qn &u)
{
    assert(n == u.n && n != 0);
    for (int i = 0; i < n; i++)
        v[i] -= u.v[i];
}

void Qn::operator=(const Qn &u)
{
    assert(n == u.n && n != 0);
    for (int i = 0; i < n; i++)
        v[i] = u.v[i];
}

void QnIP(Qn &v, int n)
{
    integer t;
    int i = 0;
    while (true)
    {
        std::cin >> t;
        v.v[i] = t;
        i++;
        if (std::cin.peek() == '\n' || i == n)
            break;
    }
}

/****************************** MATRIX **********************************/

void init(matrix &x, int n, int m)
{
    assert(x.size() == 0);
    for (int i = 0; i < m; i++)
        x.push_back(Qn(n));
}

void MatrixIp(matrix &ret)
{
    ret.clear();
    int n, m;
    std::cout << "Enter nrows and ncols: ";
    std::cin >> n >> m;
    integer x;
    Qn t(n);
    for (int i = 0; i < m; i++)
    {
        std::cout << "Enter vector " << i + 1 << " : ";
        QnIP(t, n);
        ret.push_back(t);
    }
}

void MatrixOp(matrix &ret, std::string s, int k)
{
    int n = ret[0].n, m = ret.size();
    if (s != "")
        std::cout << s << std::endl;
    integer x;
    std::cout << "[\n";
    for (int i = 0; i < n; i++)
    {
        std::cout << "[ ";
        for (int j = 0; j < m; j++)
        {
            if (k == 1)
                std::cout << ret[j].v[i] << " ";
            else
                std::cout << ret[j].v[i].numerator() << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "]\n";
}

/****************************** LLL **********************************/

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
            std::swap(basis[k], basis[k - 1]);
            GramSchmidt(basis, mu, red);
            k = std::max(k - 1, 1);
        }
    }
}
