#!/usr/bin/env bash
if g++ -std=c++0x -pthread main.cpp cxstring.cpp -o ../../simpleconv
then
    echo make_finished
fi
