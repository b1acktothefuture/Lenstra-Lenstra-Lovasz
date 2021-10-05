# LLL++
## C++ implementaion of LLL basis reduction algorithm

## Dependencies
* Boost : https://www.boost.org/
* GMP : https://gmplib.org/

## Usage

### Example 
```c++
// foo.cpp
#include "lll.h"
....
....
int main(){ 
    matrix x;
    MatrixIp(x); // takes basis for lattice as input

    lll(x, fraction(3, 4)); // delta: 𝛿 in Lovasz condition, generally 0.75
    // basis matrix is modified to store the reduced basis
    
    MatrixOp(x, "\nReduced Basis:"); // prints reduced basis to cout
}
```
### Compile
```bash
g++ foo.cpp lll.cpp -lgmp
```
#### Output
```bash
Enter nrows and ncols: 3 3
Enter Basis 1 : 1 1 1
Enter Basis 2 : -1 0 2
Enter Basis 3 : 3 5 6

Reduced Basis:
[
[ 0 1 -2 ]
[ 1 0 0 ]
[ 0 1 1 ]
]
```
## References
* https://web.eecs.umich.edu/~cpeikert/lic15/lec03.pdf
* https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm

