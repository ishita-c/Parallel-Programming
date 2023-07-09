#!/bin/bash

for (( i=0; i<=38; i++ ))
do
  input_file1="tcs/input${i}/input${i}a"
  input_file2="tcs/input${i}/input${i}b"
  output_file="tcs/input${i}/myout"
  printf "Test case ${i} \n"
  ./exec "$input_file1" "$input_file2" "$output_file"
  printf "\n"
done