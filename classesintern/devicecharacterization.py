#device characterization functions
#Juan Pablo Duarte

import numpy as np
import scipy.optimize as optimization
import supportfunctions as sp


def findvth(vg_array,ids_array,idsref):
#find vg for a given current level
  index = 0 
  diff1 = 1000
  index=0
  for ids in ids_array:
  #check1
    if abs(ids-idsref)<diff1:
      indexids1 = index
      diff1 = abs(ids-idsref)
    index+=1 
  if not ((indexids1==0) or (indexids1==(len(vg_array)-1))):
    diff2 = abs(ids_array[indexids1-1]-idsref)
    diff3 = abs(ids_array[indexids1+1]-idsref)   
  
    if(diff2<diff3):
      indexids2 = indexids1-1
    else:
      indexids2 = indexids1+1  

    vt=0.0259 
    I1=abs(ids_array[indexids1])
    I2=abs(ids_array[indexids2] )
    Vg1 = vg_array[indexids1]  
    Vg2 = vg_array[indexids2]  
    n=((Vg2-Vg1)/vt)/(np.log(I2/I1))
    Vth = n*vt*np.log(idsref/I1)+Vg1
  else:
    Vth = vg_array[indexids1] 
  
  if Vth>max(vg_array):
    Vth = max(vg_array)
  if Vth<min(vg_array):
    Vth = min(vg_array)   
   
  return Vth
  
def Ilog(vg,nvth,I0):
  return I0*np.exp(vg/nvth)

def findss(vg_array,ids_array,vgi,vthwindow,fitwindow):
#find SS with the current at given vgi wrt current*level
  #find Io for fitting curve
  index = 0 
  deltamax = 1e10
  for vg in vg_array:
    if abs(vg-vgi)<deltamax:
      indexids1 = index
      deltamax = abs(vg-vgi)
    index+=1
  Ioff = ids_array[indexids1]

  index = 0 
  for vg in vg_array:
    if abs(vg_array[0]-vg)>fitwindow: 
      deltaindex = index
      break
    index+=1
  
  index = 0
  allSS = []
  for vg in vg_array[0:-deltaindex]:
    #print 'abs(vg-vgi): '+str(vg)+'-'+str(vgi)+'='+str(abs(vg-vgi))+' vthwindow: '+str(vthwindow)
    if abs(vg-vgi)<vthwindow:
      #print "fitting"
      indexids1 = index
      #parametersfit = optimization.curve_fit(Ilog, vg_array[index:index+deltaindex], ids_array[index:index+deltaindex], (0.0256,Ioff))
      #ssaux = parametersfit[0][0]*np.log(10) #SS calculation
      Vg2 = vg_array[index]
      Vg1 =vg_array[index+1]
      I2 = abs(ids_array[index])
      I1 = abs(ids_array[index+1]) 
      ssaux=((Vg2-Vg1))/(np.log(I2/I1))*np.log(10)
      if ssaux<60e-3:
        ssaux = 1e6
      allSS.append(ssaux)
    index+=1
  if len(allSS)>0:
    SS = min(allSS)
  else:
    SS = 0
  return  SS

def findGmax(vg_array,ids_array):
  K = sp.K_generator(np.array(vg_array),1) 
  gm = K*np.array(ids_array)
  gmax = max(gm)
  return gmax

def findRon(vg_array,ids_array,vgon,vdd):
  index = 0 
  deltamax = 1e10
  for vg in vg_array:
    if abs(vg-vgon)<deltamax:
      indexids1 = index
      deltamax = abs(vg-vgon)
    index+=1
  Ion = ids_array[indexids1]
  Ron = vdd/Ion
  return Ron
  
def findIds(vg_array,ids_array,vgref,SS) : 
  index = 0 
  deltamax = 1e10
  for vg in vg_array:
    if abs(vg-vgref)<deltamax:
      indexids1 = index
      deltamax = abs(vg-vgref)
    index+=1
  if (vg_array[indexids1]!=vg_array[0] ): 
    Iref = ids_array[indexids1]
  else:
    Iref = ids_array[0]*np.exp(vgref/(SS/np.log(10)))
  return Iref
