#!/bin/bash

for i in {1..100}
do
    echo "Running iteration $i"
    python3 generate_data_lucas.py
    echo "Finished iteration $i"
done
