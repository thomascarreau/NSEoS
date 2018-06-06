#!/usr/bin/env python2

import sys
import csv
import numpy as np

NPARAMS = 18

def readPosterior():
    posterior = []
    with open(sys.argv[1]) as fposterior:
        reader = csv.reader(fposterior, delimiter=" ")
        for column in zip(*reader):
            posterior.append(column)
    fposterior.close()
    return posterior

def calcStats(posterior):
    # calculate the weights
    chi2 = posterior[-1]
    del posterior[-1]
    w = []  
    for i in range(len(posterior[0])):
        w.append(np.exp(-float(chi2[i])/2.))
    # convert to numpy array
    p = np.array(posterior).astype(np.float)
    # average
    av = []
    avm = []
    for i in range(NPARAMS):
        av.append(np.average(p[i]))
        avm.append(np.average(p[i], weights=w))
    # deviation
    sig = []
    for i in range(NPARAMS):
        sig.append(np.std(p[i]))
    # correlation matrix
    cm = np.corrcoef(p)
    #==========================================================
    statistics = open(sys.argv[2], "w")
    for i in range(NPARAMS):
        statistics.write(str(av[i]) + " " + str(sig[i]) + "\n")
    statistics.close()
    matrix = open(sys.argv[3], "w")
    for i in range(NPARAMS):
        for j in range(NPARAMS):
            matrix.write(str(cm[i][j]) + " ")
        matrix.write("\n")
    matrix.close()

def main():
    if len(sys.argv) != 4:
        sys.exit("ERROR: syntax is 'python2 stats.py posterior.in statistics.out matrix.out'")
    calcStats(readPosterior())

if __name__ == "__main__":
    main()
