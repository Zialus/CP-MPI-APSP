# cp-mpi-apsp

## How to compile

``` bash
mkdir build && cd build
cmake .. && make
```

## How to run

``` bash
mpirun -n <nprocs> -hostfile <mycluster> floyd <inputfile>
```