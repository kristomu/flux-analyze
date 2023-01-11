#!/bin/sh

g++ mfmdecoder.cc rabinkarp.cc fluxrecord.cc -lz -lsqlite3 -O3 -Wall -ggdb
