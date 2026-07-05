#!/bin/bash

N=1000
filename="grid.cmd"

for ((i=1; i<=N; i++))
do
    echo "Run $i"
    ./run.sh "$filename"
done
