# LLL++
## C++ implementaion of LLL basis reduciton algorithm

## Dependencies
* Boost : https://www.boost.org/
* GMP : https://gmplib.org/

## Usage

``` c++
#include "lll.h"
.....
matrix basis;
MatrixIp(basis); // taking basis for lattice as input

lll(basis,delta);  // delta: 𝛿 in Lovasz condition, default 0.75
// basis matrix is modified to store the reduced basis

MatrixOp(basis); // print reduced basis to cout
.....
```
### Compile
```bash
g++ foo.cpp lll.cpp -lgmp
```

## References
* https://web.eecs.umich.edu/~cpeikert/lic15/lec03.pdf
* https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm

