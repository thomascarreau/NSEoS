#!/usr/bin/env python2

import sys
import csv
import math
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

def calculateWeigths(posterior):
    chi2 = posterior[-1]
    del posterior[-1]
    weights = []  
    for i in range(len(posterior[0])):
        weights.append(np.exp(-float(chi2[i])/2.))
    return weights

def calcWeightedAverageAndDeviation(values, weights=None):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))

def calcWeightedCorrelationCoefficient(values1, values2, weights=None):
    avg_and_std1 = calcWeightedAverageAndDeviation(values1, weights)  
    avg_and_std2 = calcWeightedAverageAndDeviation(values2, weights) 
    corrcoeff = np.average((values1-avg_and_std1[0])\
            *(values2 -avg_and_std2[0]), weights=weights)\
            /avg_and_std1[1]/avg_and_std2[1]
    return corrcoeff

def calculateCorrelationMatrix(posterior, weights=None):
    m = [[], [], [], [], [], [], [], [], [], 
            [], [], [], [], [], [], [], [], []]
    for param1 in range(NPARAMS):
        for param2 in range(NPARAMS):
            m[param1].append(calcWeightedCorrelationCoefficient(posterior[param1],\
                    posterior[param2]))
    return m

def printStatistics(posterior, weights=None):
    statistics = open(sys.argv[2], "w")
    for param in range(NPARAMS):
        s = calcWeightedAverageAndDeviation(posterior[param], weights)
        statistics.write(str(s[0]) + " " + str(s[1]) + "\n")
    statistics.close()

def printCorrelationMatrix(m):
    matrix = open(sys.argv[3], "w")
    for param1 in range(NPARAMS):
        for param2 in range(NPARAMS):
            matrix.write(str(m[param1][param2]) + " ")
        matrix.write("\n")
    matrix.close()

def main():
    if len(sys.argv) != 4:
        sys.exit("ERROR: syntax is 'python2 stats.py posterior.in \
                statistics.out matrix.out'")
    posterior = readPosterior()
    w = calculateWeigths(posterior)
    p = np.array(posterior).astype(np.float)
    printStatistics(p)
    printCorrelationMatrix(calculateCorrelationMatrix(p))

if __name__ == "__main__":
    main()
