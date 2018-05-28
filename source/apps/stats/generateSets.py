#!/usr/bin/env python

from random import randint

# output file format:
# lasat0 nsat0 ksat0 qsat0 zsat0 jsym0 lsym0 ksym0 qsym0 zsym0 ms dms b p

def main():
    p = [2, 3, 4]
    for i in range(500):
        lasat0 = float(randint(150,170))/(-10.0)
        nsat0 = float(randint(150,170))/1000.0
        ksat0 = float(randint(190,270))
        qsat0 = float(randint(-1000,1000))
        zsat0 = float(randint(-3000,3000))
        jsym0 = float(randint(260,380))/10.0
        lsym0 = float(randint(10,80))
        ksym0 = float(randint(-400,200))
        qsym0 = float(randint(-2000,2000))
        zsym0 = float(randint(-5000,5000))
        ms = float(randint(60,80))/100.0
        dms = float(randint(0,20))/100.0
        b = float(randint(10,100))/10.0
        for pval in p:
            print lasat0, nsat0, ksat0, qsat0, zsat0,  \
                    jsym0, lsym0, ksym0, qsym0, zsym0, \
                    ms, dms, b, pval

if __name__ == "__main__":
    main()
