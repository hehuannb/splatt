#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "linux" && "$CC" == "gcc" ]]; then
  gcov `find build -name '*.gcno'`
fi

