#!/usr/bin/env bash

python2 randomSets.py > ../../input/rng_sets.data
make
./posterior ../../input/jm_sets.data posterior_ld.out
./posterior ../../input/mw_sets.data posterior_hd.out
./posterior ../../input/rng_sets.data posterior_rng.out
