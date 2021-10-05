#include "src/lll.h"

void test()
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
    MatrixOp(x, "Reduced Basis:");
}

void user_test()
{
    matrix x;
    MatrixIp(x);
    lll(x, fraction(3, 4));
    MatrixOp(x, "\nReduced Basis:");
}

int main()
{
    user_test();
}