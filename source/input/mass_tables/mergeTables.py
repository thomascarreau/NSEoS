import sys
import numpy as np

def CheckIfInTables(table1, table2, zz, nn):
    '''Return mass excess for zz protons and nn neutrons of table1, or table2 
    if not found in table1. Return None if not found in both tables.'''
    check = 0
    for i in range(len(table1[0])):
        if table1[0][i] == zz and table1[1][i] == nn:
            check = 1
            return table1[2][i]
    for i in range(len(table2[0])):
        if table2[0][i] == zz and table2[1][i] == nn:
            return table2[2][i]

def MergeTables(path_of_table1, path_of_table2, path_of_merged_table):
    '''Merge the two input mass tables. Mass excess of table1 is kept in 
    priority. Reference table for zz and nn is HFB-24.'''
    ref = np.loadtxt('hfb24.data', unpack=True) # reference table
    table1 = np.loadtxt(path_of_table1, unpack=True)
    table2 = np.loadtxt(path_of_table2, unpack=True)
    hfb24 = np.loadtxt('hfb24.data', unpack=True)
    merged_table = open(path_of_merged_table, 'w')
    for i in range(len(ref[0])):
        deps = CheckIfInTables(table1, table2, ref[0][i], ref[1][i])
        if deps != None:
            merged_table.write(str(int(ref[0][i])) + ' ' + str(int(ref[1][i])) 
                    + ' ' + str(deps) + '\n')
    merged_table.close()

def main():
    path_of_table1 = sys.argv[1]
    path_of_table2 = sys.argv[2]
    path_of_merged_table = sys.argv[3]
    MergeTables(path_of_table1, path_of_table2, path_of_merged_table)

if __name__=='__main__':
    main()
