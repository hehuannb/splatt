#!/bin/bash

gcov --version

for f in $(find build -name '*.gcno');
do
  gcov ${f}
done
#gcov -o $(find build -name '*.gcno')

