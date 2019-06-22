# PLOT_EIG: Plot Eigenvalue Estimates from TMP
#
# Author: Thomas R. Cameron
# Date: 6/22/2019
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# open data file
f = open('../data_files/eigvals.csv')
lineList = f.readlines()
itnum = len(lineList)
eignum = len(lineList[0].split(","))
eigvals = [np.zeros(eignum,dtype=complex) for k in range(itnum)]
# store eigenvalue estimates
for k in range(itnum):
    row = lineList.pop(0)
    row = row.split(",")
    for l in range(eignum):
        eigvals[k][l] = eval(row[l])
# figure setup
colors = mpl.cm.Blues(np.linspace(0.25,1,itnum))
# plot eigenvalues
plt.scatter(np.real(eigvals[0]),np.imag(eigvals[0]),color=colors[0],label='InitEst')
#for k in range(1,itnum-1):
#    plt.scatter(np.real(eigvals[k]),np.imag(eigvals[k]),color=colors[k])
plt.scatter(np.real(eigvals[itnum-1]),np.imag(eigvals[itnum-1]),color=colors[itnum-1],label='FinEst')
plt.xlabel('Real')
plt.ylabel('Imaginary')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=2, mode="expand", borderaxespad=0.)
plt.show()