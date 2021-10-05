#include "../src/lll.h"

/* 
checks if given instance has low density. 
density = n / log2 ( max (arr) )
density condition max(arr) >= 2**(n*n/2 + epsilon), hence acceptable value < 2/n 
*/



bool checkDensity(integer n, integer X){
    fraction t(n*n,2);
    fraction x(0);
    for(;X!=0;X>>=1) x++;
    return x>t;
}

integer checkSolution(matrix& reducedBasis,Zn arr,integer sum,bool flip = 0){
    integer t(0),ret(0);

    for(int i = 0;i<reducedBasis.size();i++){
        t = 0;
        ret = 0;
        for(int j = 0;j<arr.size();j++){
            if(reducedBasis[i].v[j]!=integer(0)){
                t += arr[j];
                if(!flip) ret += integer(1)<<j;
            }
            else if(flip) ret += integer(1)<<j;
        }
        if(t == sum) return ret;
    }
    return -1;
}


integer lowDensityKnapsack(Zn arr,integer sum,bool flip = 0){

    // if sum is zero return 0

    if(sum == 0)return 0; 

    int n = arr.size();
    integer total(0),largestX(0);
    for(int i = 0;i<n;i++){
        total += arr[i]; // sum(arr)
        if(arr[i] > largestX) largestX = arr[i]; // max(arr)
    }

    // check density condition and solvability  

    if(sum >total){
        std::cout << "Invalid Instance: given S is greater than sum of array";
        return -1;
    }

    if(!checkDensity(arr.size(),largestX)) {
        std::cout<<"Invalid Instance: high density"; 
        return -1; 
    }
    
    // Solves dubset dum for SUM(arr)-sumand returns negtion fof result
    if(2*sum< total) return (lowDensityKnapsack(arr,total-sum,1));

/*
    The lattice basis matrix (n+1 x n+1) is:
    basis = [I   O] 
            [b*A S]

    I = nxn identity matrix
    A = 1xn input arr
    O = nx1 zero vector
    S = given sum
    b = ceil(sqrt(n x 2^n))
*/

    matrix basis;
    init(basis,n+1,n+1);

    integer bSquare(n*(integer(1)<<n));
    floating y = boost::multiprecision::sqrt(bSquare);
    integer b(y+1);
    
    for(int i =0;i<n;i++){
        basis[i].v[i] = 1;
        basis[i].v[n] = -1*b*arr[i];
    }
    basis[n].v[n] = sum;
    
    // Basis Reduction using LLL, delta = 0.75
    // std::cout << "Reducing Lattice Basis...\n";
    lll(basis,fraction(3,4)); 

    return checkSolution(basis,arr,sum,flip);
}




int main(){
    std::cout << "Enter n (number of elements): ";
    int x;
    std::cin>> x;
    Zn v;
    for(int i = 0;i<x;i++){
        integer t;
        std::cin>> t;
        v.push_back(t);
    }

    std::cout<<lowDensityKnapsack(v,v[1]+v[2]+v[3])<<std::endl;
}