#!/bin/bash

num=$1

rm pico_prod_random.list
rm runNumber_prod_random.list

shuf -n ${num} pico_prod_sorted.list > pico_prod_random.list
awk -F/ '{print $(NF-1), $0}' pico_prod_random.list | cut -d ' ' -f 1 | sort | uniq > runNumber_prod_random.list
