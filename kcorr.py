from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.io import ascii
import decimal
from scipy import constants
from pylab import *
import pylab as py

#Global 

redshiftsize=2.0
redshiftbins=20
maxy=40

#Array containing all redshifts bins/redshifts
 
dz=float(redshiftsize/redshiftbins)
bins=arange(0, redshiftsize+dz, dz)
print(bins, bins.shape[0], bins[3])

#Matrix Containing Partial SUMS

FINAL=np.zeros((bins.shape[0],4))

#Open SED: SED will expand 
file1=ascii.read('m51_int')

#Wavelength and SED Flux columns
asedwave, aflux=file1['wave'], file1['flux']
sedwave, flux=np.array(asedwave), np.array(aflux)

#Want Wavelengths to have exactly two decimal places --> round
#Number of rows in SED file
sedrows=sedwave.shape[0]

#Create reference SED Matrix (Biggest One)
SED=np.zeros((sedrows,2))

#Modifiable SED matrix
zSED=np.zeros((sedrows,2))

#Fill Matrix with rounded wavelength values
#First column 1 and then column 2 

count=0
while (count<sedrows):
    SED[count][0]=round(float(sedwave[count]),2)
    zSED[count][0]=round(float(sedwave[count]),2)
    count=count+1
 
count=0
while (count<sedrows):
    SED[count][1]=flux[count]
    zSED[count][1]=flux[count]
    count=count+1

#Open FILTER file: the filter doesn't change
file2=ascii.read('24_int')

#Wavelength and Transmission columns
aoriginalwave, atrans=file2['wave'], file2['final']
originalwave, trans=np.array(aoriginalwave), np.array(atrans)

#Number of rows
rows=originalwave.shape[0]
 
#Lowest filter value
FILTERlowest=originalwave[0]

######SEDS: REDSHIFT STRETCHING######

####REDSHIFT COUNTER###########
for z in range(bins.shape[0]):#
###############################
  
  #######################################
  #SAVE ASCII FILE WITH STRETCHED VALUES#
  #######################################
  
  stretchedfile=open('stretchz='+str(bins[z]), 'w')
  stretchedfile.write('#zwave flux \n\n')
  
  #CREATE STRETCHED SED VALUES
  for i in range(sedrows):
    zSED[i][0]=round(float(SED[i][0]*(1+bins[z])),2)
  
  for i in range(sedrows):
    result="{0} {1}".format(zSED[i][0], zSED[i][1])
    stretchedfile.write(result+'\n')
  
  stretchedfile.close()
  
  #LOWEST/HIGHEST VALUES OF zSED FILE  
  zSEDlow=round(float(zSED[0][0]+0.01),2)
  zSEDhigh=round(float(zSED[sedrows-1][0]),2)
  #print(zSED, zSEDlow, zSEDhigh) DEBUG
  
  ##################
  #INTERPOLATE zSED#
  ##################
  #############################################
  #READ SAVED ASCII FILE WITH STRETCHED VALUES#
  #############################################
  
  file3=ascii.read('stretchz='+str(bins[z]))
  aswave, asflux=file3['zwave'], file3['flux']
  swave, sflux=np.array(aswave), np.array(asflux)
  #print(aswave, asflux)
  
  inter=interp1d(swave, sflux, kind='linear')
  
  #######################################
  #SAVE INTERPOLATED zSED TO ASCII FILES#
  #######################################
  
  #Minimum and Maximum values for interpolation
  intercount=round(float(zSED[0][0]+0.01),2)
  arraymax=round(float(zSED[sedrows-1][0]),2)
  
  results1=open('interz='+str(bins[z]),'w')
  results1.write("#zwave flux \n\n")
  
  while (intercount < arraymax):
    interpolated=inter(intercount)
    result="{0} {1}".format(intercount, inter(intercount))
    results1.write(result+'\n')
    intercount=intercount+0.01
    
  results1.close()
  
  ###### VALUE MATCHING (FILTER -> SED)######
 
  ##################################################
  #NUMBER OF ROWS BETWEEN zSEDlow and FILTERlowest #
  ##################################################
  
  file4=ascii.read('interz='+str(bins[z]))
  aintwave, aintflux=file4['zwave'], file4['flux']
  intwave, intflux=np.array(aintwave), np.array(aintflux)
  
  diff=round(float((FILTERlowest-intwave[0])*100),2)
  #print(diff, FILTERlowest, intwave[0], intwave[diff]) DEBUG
  
  ##############################################################
  #FILTERlowest value is found at SED[diff][0]: Multiply Values#
  ##############################################################
  
  SUM=0
  for width in range(rows):
      SUM=SUM+(intflux[diff+width]*trans[width])
  
  print('z='+str(bins[z]), 'SUM='+str(SUM))
  
  #############################
  #STORE VALUES IN MATRIX SUMS#
  #############################
  
  ####################
  #FIRST K VALUE IS 1#
  ####################
  
  if z == 0:
    FINAL[z][0]=bins[z]
    FINAL[z][1]=SUM
    FINAL[z][2]=1
    FINAL[z][3]=1
  else:
    FINAL[z][0]=bins[z]
    FINAL[z][1]=SUM
    FINAL[z][2]=FINAL[0][1]/FINAL[z][1]
    FINAL[z][3]=FINAL[z][2]/(1+bins[z])

print(FINAL)

#####################################################   
#ASCII FILE TO STORE FINAL DATA POINTS (file: FINAL)#
#####################################################

final=open('FINAL','w')
final.write("#z SUM K K/(1+z) \n\n")

for z in range(bins.shape[0]):
 result="{0} {1} {2} {3}".format(FINAL[z][0], FINAL[z][1], FINAL[z][2], FINAL[z][3])
 final.write(result+'\n')
    
final.close()

#############################
#PLOT DATA POINTS FROM FINAL#
#############################

#############################
#READ DATA POINTS FROM FINAL#
#############################

file5=ascii.read('FINAL')
aredshift, acorr=file5['z'], file5['K/(1+z)']
redshift, corr=np.array(aredshift), np.array(acorr)

####################################
#INTERPOLATE DATA POINTS FROM FINAL#
####################################

interfinal=interp1d(redshift, corr, kind='cubic')

######################################
#SAVE INTERPOLATED DATA IN interFINAL#
######################################

intfinal=open('interFINAL','w')
intfinal.write("#z K/(1+z) \n\n")

z=0
while (z <= redshiftsize):
  interpolated=interfinal(z)
  result="{0} {1}".format(z, interfinal(z))
  intfinal.write(result+'\n')
  z=z+0.01
    
intfinal.close()

##############
#PLOT RESULTS#
##############

plt.plot(redshift, corr, 'o')
plotx=np.linspace(0, 2, num=100, endpoint=True)
plt.plot(plotx, interfinal(plotx), '--')
plt.yscale('log')
py.ylim([0,maxy])
#plt.yticks(np.arange(min(interfinal(plotx)), max(interfinal(plotx))+1, 5.0))
py.ylabel("K/(1+z)")
py.xlabel("z")
plt.show()
