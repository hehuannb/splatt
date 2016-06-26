#!/bin/bash

find build -name '*.gcno'
gcov `find build -name '*.gcno'`

