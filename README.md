# cp-mpi-apsp

[language-badge]: https://img.shields.io/badge/C%2B%2B-11-red.svg
[platform-badge]: https://img.shields.io/badge/Platforms-OS%20X%20--%20Linux-green.svg?style=flat

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c8e2489cd4f849bcbb60a33049f2ab02)](https://www.codacy.com/app/raulmendesferreira/cp-mpi-apsp?utm_source=github.com&utm_medium=referral&utm_content=Zialus/CP-MPI-APSP&utm_campaign=badger)
[![Build Status](https://travis-ci.org/Zialus/CP-MPI-APSP.svg?branch=master)](https://travis-ci.org/Zialus/CP-MPI-APSP)
[![Build status](https://ci.appveyor.com/api/projects/status/25scbj2qhkisyytn?svg=true)](https://ci.appveyor.com/project/Zialus/cp-mpi-apsp)

## How to compile

``` bash
mkdir build && cd build
cmake .. && make
```

## How to run

``` bash
mpirun -n <nprocs> -hostfile <mycluster> floyd <inputfile> <outputfile>
```
