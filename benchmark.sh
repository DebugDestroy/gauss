#!/bin/bash

N=100
filename="grid.cmd"

for ((i=1; i<=N; i++))
do
    echo "Run $i"
    ./run.sh "$filename"
done
