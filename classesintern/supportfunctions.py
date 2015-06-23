#aux functions
import numpy as np
from numpy.matlib import repmat
from scipy.misc import factorial
from scipy import sparse
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm

def meshgrid2(*arrs):
#this generate array similar to ngrid in matlab for many input vectors
    arrs = tuple(reversed(arrs))  #edit
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []    
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = np.asarray(arr).reshape(slc) #array Convert the input to an array
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j) 
        ans.append(arr2)
        
    #print "type ans: ", type(ans) #ans is a list
    i=0
    totallength = np.product(ans[0].shape)
    for variablearray in ans[::-1]:
        ans[i] = np.concatenate(variablearray.reshape(totallength,1))
        i=i+1
    return ans
    
def inplace_change(filename, old_string, new_string):
        s=open(filename).read()
        if old_string in s:
                #print 'Changing "{old_string}" to "{new_string}"'.format(**locals())
                s=s.replace(old_string, new_string)
                f=open(filename, 'w')
                f.write(s)
                f.flush()
                f.close()
        else:
                print 'No occurances of "{old_string}" found.'.format(**locals())
                

#Derivative Support functions
def mkfdstencil(x,xbar,k):
#this funtion is sue to create finite diference method matrix
  maxorder            = len(x)
  h_matrix            = repmat(np.transpose(x)-xbar,maxorder,1)
  powerfactor_matrix  = np.transpose(repmat(np.arange(0,maxorder),maxorder,1))
  factorialindex      = np.transpose(repmat(factorial(np.arange(0,maxorder)),maxorder,1))
  taylormatrix        = h_matrix ** powerfactor_matrix /factorialindex
  derivativeindex     = np.zeros(maxorder)
  derivativeindex[k]  = 1 
  u = np.linalg.solve(taylormatrix,derivativeindex)
  return u

def K_generator(x,order):
#this return matrix to find the derivative, x is the variable to be derived and order is the derivative order
  N=len(x);
  K = lil_matrix((N, N))
  K[0,:6]=mkfdstencil(x[0:6],x[0],order)
  K[1,:6]=mkfdstencil(x[0:6],x[1],order)
  K[2,:6]=mkfdstencil(x[0:6],x[2],order)

  i=3
  for xbar in x[3:-3]:
    #print i
    K[i,i-3:i+3]=mkfdstencil(x[i-3:i+3],xbar,order)
    i+=1
  #print i
  K[i,-7:-1]=mkfdstencil(x[-7:-1],x[-3],order)  
  i+=1
  K[i,-7:-1]=mkfdstencil(x[-7:-1],x[-2],order)
  i+=1  
  K[i,-7:-1]=mkfdstencil(x[-7:-1],x[-1],order)  
  return K.tocsr()

#########################################################################
#arrange X and Y matrices
def rearrangearray(arrayXa,elementpercylce,numberelement):
#this function reshpae array to be printing 
  arrayXb = arrayXa.reshape((elementpercylce, len(arrayXa)/elementpercylce))#arrayXa.reshape(( len(arrayXa)/elementpercylce,elementpercylce))#
  arrayXc = np.transpose(arrayXb)
  arrayXd = arrayXc.reshape((len(arrayXa)/numberelement,numberelement))#arrayXc.reshape((numberelement,len(arrayXa)/numberelement))#
  arrayXe = np.transpose(arrayXd)
  return arrayXe

def findelementpercylce(arrayaux):
#this function return the number of entire sequence in the initial array
  #first find number of elements repeated next to each other
  firstelement = arrayaux[0]
  lengthaux = len(arrayaux)
  flag=1
  i=1

  elementpercylce = 0
  while ( (i<(lengthaux+1)) and flag):
    elementpercylce = i-1

    if (abs(arrayaux[i]-firstelement)>0): #%TODO: check abs condition
      flag=0
    i=i+1

  elementpercylce = elementpercylce+1
  
  #this return number of time a sequence repeat
  indexes = []
  b = arrayaux[0:elementpercylce]
  for i in range(len(arrayaux)-len(b)+1):
    if sum(abs(arrayaux[i:i+len(b)] - b))<1e-15:
      indexes.append((i, i+len(b)))
  return len(indexes)
  #return elementpercylce                
