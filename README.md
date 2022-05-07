# Parallelized Genetic Algorithm 
Authors : [Mohit Gayatri](https://github.com/mgayatri77) and [Sridhar Pandian Arunachalam](https://github.com/SridharPandian). 

This repository contains the code for our High Performance Computing Project for Spring 2022. This is a C/C++ implementation of a standard Genetic Algorithm with various objective functions. This was parallelized using OpenMP (CPU based parallelization).

## Prerequisites
- C/C++ version 11 or greater
- OpenMP version 4.5

## Running the algorithm
Use the following command to compile and run a serial version of the islands algorithm:
```
g++ -O3 -std=c++11 main.cpp -o main_serial && ./main_serial
```
Use the following command to compile and run a parallelized version of the islands algorithm:
```
g++ -O3 -std=c++11 main.cpp -o main_parallel && ./main_parallel
```