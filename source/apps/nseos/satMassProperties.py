#!/usr/bin/env python2

from math import sqrt

def CalcKs(m):
    return 1./m - 1.

def CalcKv(m, dm):
    ks = CalcKs(m)
    if (dm == 0.0):
        return (sqrt(dm*dm*(1.+ks)*(1.+ks) + 1.) + dm*ks - 1.)/dm
    else:
        return ks - 0.5*dm*(1.+ks)*(1.+ks)

def CalcBardel(m, dm):
    ks = CalcKs(m)
    kv = CalcKv(m, dm)
    return ks - kv

def DoComparisonWithPaper():
    mass_satdata = [[0.67,0.26,"ratp   ",0.78], [0.80,0.03,"bsk14  ",0.28], 
            [0.80,0.04,"bsk16  ",0.28], [0.80,0.04,"bsk17  ",0.28], 
            [0.69,-0.19,"sly4   ",0.25], [0.69,0.21,"nrapr  ",0.66], 
            [0.70,-0.47,"sly230a",0.00], [0.69,-0.19,"sly230b",0.25], 
            [0.90,0.09,"sko    ",0.17], [0.66,-0.06,"ddme1  ",0.45], 
            [0.65,-0.06,"ddme2  ",0.47], [0.67,-0.08,"nl3    ",0.40],
            [0.71,-0.09,"tm1    ",0.32]]
    print "set, my kv, kv of paper"
    print "======================="
    for i in range(len(mass_satdata)):
        kv = CalcKv(mass_satdata[i][0], mass_satdata[i][1])
        print mass_satdata[i][2], " ", round(kv,2), " ", mass_satdata[i][3]
    print ""

def main():
    DoComparisonWithPaper()
    keeplooping = 1
    while (keeplooping == 1):
        m = float(raw_input("m = "))
        dm = float(raw_input("dm = "))
        print "ks     =", CalcKs(m)
        print "kv     =", CalcKv(m, dm)
        print "bardel =", CalcBardel(m, dm)
        doloop = raw_input("Enter y or Y to continue... ")
        if (doloop not in ["y", "Y"]):
            keeplooping = 0
        print ""

if __name__ == "__main__":
    main()
