#!/bin/bash

N=10
filename="rrt.cmd"

for ((i=1; i<=N; i++))
do
    echo "Run $i"
    ./run.sh "$filename"
done
