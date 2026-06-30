#!/bin/bash

N=1000
filename="project1.cmd"

for ((i=1; i<=N; i++))
do
    echo "Run $i"
    ./run.sh "$filename"
done
