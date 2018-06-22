#!/usr/bin/env python2

from random import randint

# output file format:
# lasat0 nsat0 ksat0 qsat0 zsat0 jsym0 lsym0 ksym0 qsym0 zsym0 ms dms b

def GenerateSet():
    ''' Randomly generate a set of empirical parameters
    given a uniform distribution '''
    print \
            float(randint(150,170))/(-10.0), \
            float(randint(150,170))/1000.0, \
            float(randint(190,270)), \
            float(randint(-1000,1000)), \
            float(randint(-3000,3000)), \
            float(randint(260,380))/10.0, \
            float(randint(10,80)), \
            float(randint(-400,200)), \
            float(randint(-2000,2000)), \
            float(randint(-5000,5000)), \
            float(randint(60,80))/100.0, \
            float(randint(0,20))/100.0, \
            float(randint(10,100))/10.0

def main():
    for i in range(1000):
        GenerateSet()

if __name__ == "__main__":
    main()
