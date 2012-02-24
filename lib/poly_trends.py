#!/frodo/shared/epd/bin/python
import numpy as np
import sys
import os

#constants:
nTRs = int(sys.argv[1])
nPoly = int(sys.argv[2])


## Create polynomial regressors centered on 0 and bounded by -1 and 1
x = np.arange(nTRs)
num_pol = range(nPoly)
y = np.ones((len(num_pol),len(x)))


for i in num_pol:
    y[i,:] = (x - (np.max(x)/2)) **(i+1)
    y[i,:] = y[i,:] - np.mean(y[i,:])
    y[i,:] = y[i,:]/np.max(y[i,:])    




## Print out text file for each polynomial to be used as a regressor
for i in num_pol:
    np.savetxt('poly_detrend_' + str(i+1) + '.txt',y[i],fmt='%0.02f')
