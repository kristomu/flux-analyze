#!/bin/sh

g++ mfmdecoder.cc rabinkarp.cc -lz -lsqlite3 -O3 -Wall -ggdb
