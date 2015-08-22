from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.io import ascii
import decimal

#Global 

centre=100
resolution=0.25
redshiftsize=3.0
redshiftbins=12

#Open SED

file1=ascii.read('m51.sed')
#print(file1)
colthree=raw_input('col1\n')
colfour=raw_input('col2\n')

#Wavelength and SED Flux columns
asedwave, aflux=file1[colthree], file1[colfour]
sedwave, flux=np.array(asedwave), np.array(aflux)

#Want Wavelengths to have exactly two decimal places --> round
#Number of rows in SED file
sedrows=sedwave.shape[0]

#Create SED Matrix (Biggest One)

SED=np.zeros((sedrows,2))

#Fill Matrix with rounded wavelength values
#First column 1 and then column 2 

count=0
while (count<sedrows):
    SED[count][0]=round(float(sedwave[count]),2)
    count=count+1

count=0

while (count<sedrows):
    SED[count][1]=flux[count]
    count=count+1

#print(SED[0][0])
#print(SED.shape)
#Want to create a matrix with all shifted wavelengths
#Open FILTER file

file2=ascii.read('herschel100.filter')
#print(file2)
colone=raw_input('col1\n')
coltwo=raw_input('col2\n')

#Wavelength and Transmission columns
aoriginalwave, atrans=file2[colone], file2[coltwo]
originalwave, trans=np.array(aoriginalwave), np.array(atrans)

#Number of rows
rows=originalwave.shape[0]

#Lowest filter value
lowest=originalwave[0]

#Array containing all redshifts bins/redshifts

dz=float(3.0/redshiftbins)
print(dz)
bins=arange(0, redshiftsize+dz, dz)
print(bins)

FILTER=np.zeros((rows,2))

#Array containing SED*TRANSMISSION results

RESULT=np.zeros((rows,4))

#Get maximum and minimum values for all redshifts

#----------------- LOOP TO COMPUTE ALL REDSHIFTS (ADD LATER) -----------------

#i can go from 0 to redshiftbins
        
for z in range(redshiftbins):
  #Central Wavelength change
  newcentre=centre/(1+bins[z])
  print(bins[z], newcentre)

  #newcentre=centre/(1+3.5)
  diff=centre-newcentre

  #Lowest filter's value shifts 
  newlowest=round(float(lowest-diff),2)

  #print(newcentre, diff, newlowest)
  #Create a matrix that contains redshifted wavelengths, but with same Transmission
  #First: Modified Wavelengths

  count2=0

  while (count2<rows):
    FILTER[count2][0]=newlowest+count2*0.01
    count2=count2+1
    
  #Second: Trasmission Function

  count2=0

  while (count2<rows):
    FILTER[count2][1]=trans[count2]
    count2=count2+1  

  results=open('redshiftedfilter', 'w')
  results.write("#wave trans\n\n") 

  for i in range(rows):
    result="{0}, {1}".format(FILTER[i][0], FILTER[i][1])
    results.write(result+'\n')

  #Find the first value from the FILTER matrix inside the SED matrix
  if FILTER[0][0] < 0: #When values in the filter file are negative we have a problem
    negative=1		
    while negative < rows: #Skip first value and go to the next one until find a positive one
      if FILTER[negative][0] < 0:
	#print(negative, FILTER[negative][0])
	negative=negative+1
	#print(negative, FILTER[negative][0])
      else:
	break
  
    ncount=negative

    while ncount < rows:
      for first in range(sedrows):
	fsed=round(float(SED[first][0]),2)
	ffilter=round(float(FILTER[ncount][0]),2)
	zero=round(float(fsed-ffilter),2)
	#print(fsed, ffilter, zero, ncount)
	if zero == 0:
	  #print(fsed,ffilter, zero, ncount)
	  break
      if zero == 0:
	#print(fsed,ffilter, zero, ncount)
	break
      else:
	ncount=ncount+1
    #print(ncount, FILTER[ncount][0])
      
  else:
    for first in range(sedrows):
      fsed=round(float(SED[first][0]),2)
      ffilter=round(float(FILTER[0][0]),2)
      zero=round(float(fsed-ffilter),2)
      #print(fsed, ffilter, zero, first)
      if zero == 0:
	#print(fsed, ffilter, zero, first)
	break


  common=0 #counter for both RESULT and FILTER (if the first value in the filter file is not negative)
  fsum=0 #total flux (sum of all TRANSMISSION*FILTER values)


  #ASCII FILE
  results=open('results.txt', 'w')
  results.write("#wavesed trans flux total\n\n") 

  if FILTER[0][0] < 0: #Skip the first n negative values
		       #IMPORTANT: when first value of FILTER matrix is negative, then we know SED[ncount][0]=FILTER[0][0]
  
    negativecounter=ncount
    common=0
    while negativecounter < rows:
  
      RESULT[common][0]=SED[common][0] #wavelength column from SED
      RESULT[common][1]=FILTER[negativecounter][1] #transmission value
      RESULT[common][2]=SED[common][1] #SED value
      RESULT[common][3]=FILTER[negativecounter][1]*SED[common][1] #SED*TRANSMISSION
      result="{0}, {1}, {2}, {3}".format(RESULT[common][0], RESULT[common][1], RESULT[common][2], RESULT[common][3])
      results.write(result+'\n')
      fsum=fsum+RESULT[common][3]
      negativecounter=negativecounter+1
      common=common+1

    results.close()
    print('z=', bins[z], 'total flux', fsum)

  else:
  
    cfirst=first #first used as a counter
    #Add "rows" to "first" to find all the FILTER values inside the SED matrix
    positivetotal=first+rows
    #print(first, rows, positivetotal)
    while cfirst < positivetotal:
  
      RESULT[common][0]=SED[cfirst][0] #wavelength column from SED
      RESULT[common][1]=FILTER[common][1] #transmission value
      RESULT[common][2]=SED[cfirst][1] #SED value
      RESULT[common][3]=FILTER[common][1]*SED[cfirst][1] #SED*TRANSMISSION
      result="{0}, {1}, {2}, {3}".format(RESULT[common][0], RESULT[common][1], RESULT[common][2], RESULT[common][3])
      results.write(result+'\n')
      fsum=fsum+RESULT[common][3]
      cfirst=cfirst+1
      common=common+1

    results.close()
    print('z=', bins[z], 'total flux', fsum)
