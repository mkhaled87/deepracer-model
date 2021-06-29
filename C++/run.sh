#!/bin/sh

g++ -c deepracer.cpp
g++ -o deepracer.exe deepracer.o
rm deepracer.o

./deepracer.exe

