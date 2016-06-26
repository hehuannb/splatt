#!/bin/bash

gcov --version
gcov $(find build -name '*.gcno')

