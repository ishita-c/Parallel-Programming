#!/bin/bash

for (( i=0; i<=38; i++ ))
do
  gold_out="tcs/input${i}/output${i}"
  my_out="tcs/input${i}/myout"
  printf "Test case ${i} \n"
  python3 checker.py -f "$my_out" "$gold_out" -e 4
  printf "\n"
done