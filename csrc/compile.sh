#!/bin/sh

g++ mfm_decoder.cc rabin_karp.cc flux_record.cc pulse_train.cc -lz -lsqlite3 -O3 -Wall -ggdb
