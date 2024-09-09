# N2HE-cosine-similarity-batch
This is a C++ open-source library which implements a FHE-based batched cosine similarity evaluation algorithm. 

We use the BFV scheme as the underlying FHE scheme and design a batched cosine similarity evaluation algorithm. The input is an encrypted n-dimension vector Enc($\vec{a}$) and $N$ encrypted vectors Enc($\vec{b}_1,\vec{b}_2,...,\vec{b}_N$).
All these vectors should have $l_2$-norm $1$. 

The algorithm will output all the $N$ cosine similarities of $\vec{a}$ and $\vec{b}_i$, $\forall i$. 


## Prerequisites
- [OpenSSL](https://www.openssl.org/)
- [hexl](https://github.com/intel/hexl)
- [SEAL](https://github.com/microsoft/SEAL)
- [Openmp](https://www.openmp.org)

## Installation
Installation on Linux:  

```
mkdir build && cd build
cmake ..
make
```

Run the test: 
```
./test
```

## License
This software is distributed under the BSD-3-Clause-Clear license. 
