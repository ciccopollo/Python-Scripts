from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.io import ascii
import decimal
import pylab
 
#Field Area

area=9000000
print(area)
s3sigma=2716.0
#Search Radius
rs=7.5


#Find n3sigma

n3sigma=s3sigma/area
print(n3sigma)

#Open catalog

#catalog with only fluxes
catalog=ascii.read('cosmos_flux_only')
afluxes=catalog['radio_flux(Jy)']
fluxes=np.array(afluxes)
#total number of rows in catalog
totrows=fluxes.shape[0]
print(totrows)

#Open matched sources file
matches=ascii.read('cosmos_matches')
aflux, ar, asubname, aradioname=matches['radio_flux(Jy)'], matches['r'], matches['Sub_id'], matches['radio_name']
flux, r, subname, radioname=np.array(aflux), np.array(ar), np.array(asubname), np.array(aradioname)
#total number of rows in catalog
mrows=flux.shape[0]

#results matrix
results=np.zeros((mrows,3))
print(mrows)

#ASCII FILE

#Ascii Header
results=open('p-values', 'w')
results.write("#Submm_Name Radio_name Flux(Jy) radius ns p\n\n") 

#For each value of the 'matched sources' file select a flux value
x=0
for x in range(mrows):
  y=0
  #print(x,y)
  while (y < totrows):
      #print(flux[x], fluxes[y])
      #print(flux[x]-fluxes[y])
      if flux[x]-fluxes[y]>=0:
	y=y+1
      else:
	lowrow=y
	sgreater=totrows-lowrow
	print("lowrow:", lowrow, "row:", x, "higher", sgreater)
	ns=float(sgreater)/area
	#print(sgreater, ns, r[x])
	ps=1-exp(-(np.pi)*ns*pow(r[x],2))
	pc=np.pi*n3sigma*pow(rs,2)
	p=1-exp(-ps*(1+log(pc/ps)))
	#rounded values
	rps=round(float(ps),5)
	rpc=round(float(pc),5)
	rp=round(float(p),5)
	result="{0} {1} {2} {3} {4} {5}".format(subname[x], radioname[x], flux[x] ,r[x], ns,rp)
	results.write(result+'\n')
	#print(rps, rpc, rp)
	break

results.close()





