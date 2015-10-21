from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.io import ascii
import decimal
import pylab

#Open catalog

#Open matched sources file
matches=ascii.read('p-values')
aflux, ar, asubname, aradioname, ap=matches['Flux(Jy)'], matches['Separation'], matches['Submm_name'], matches['Radio_name'], matches['p']
flux, r, subname, radioname, p=np.array(aflux), np.array(ar), np.array(asubname), np.array(aradioname), np.array(ap)
#total number of rows in catalog
mrows=flux.shape[0]

#ASCII FILE

#Ascii Header
results=open('robustness', 'w')
results.write("#Submm_Name Radio_name Flux(Jy) radius ns p\n\n") 

#For each value of the 'matched sources' file select a flux value
x=0
for x in range(mrows):
  if p[x]<=0.05:
    result="{0} {1} {2} {3} {4} {5}".format(subname[x], radioname[x], flux[x] ,r[x], p[x], "robust")
    results.write(result+'\n')
  elif 0.05 < p[x] < 0.2:	
    result="{0} {1} {2} {3} {4} {5}".format(subname[x], radioname[x], flux[x] ,r[x], p[x], "tentative")
    results.write(result+'\n')
  else:
    result="{0} {1} {2} {3} {4} {5}".format(subname[x], radioname[x], flux[x] ,r[x], p[x], "n/a")
    results.write(result+'\n')

results.close()





 
