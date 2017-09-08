#!/bin/bash


g++ *.cpp -Wall -Wextra -O3 -pthread -std=c++14 -o net_open_finder -static-libstdc++ -pthread -fopenmp

