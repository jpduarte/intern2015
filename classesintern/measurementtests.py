#this class is to handle data obtained from measurements
#Juan Pablo Duarte

import os
import supportfunctions as sf
import matplotlib.pyplot as plt
from numpy import loadtxt
import numpy as np

from numpy.matlib import repmat
from scipy.misc import factorial
from scipy import sparse
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from scipy import interpolate

#########################################################################
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

#########################################################################
#class definition
class mt:
  def __init__(self):#, model):
    self.version = 'v1'
    #self.model = model
    
    #defaul parameters
    self.symbol  = 'o'
    self.color = ''
    self.markerfacecolor = (1, 1, 1, 1)
    self.lw=1
    self.ylogflag = 0
    self.xlogflag = 0
    self.derivativeorder = 0
    self.markersize= 10
    self.filetype = 'SAS'
    self.devices = {}
    self.device_template = {}
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    if type(value) == type(''):
      exec "self."+name+' = '+'\''+value+'\''
      #print "self."+name+' = '+'\''+value+'\''
    else:
      exec "self."+name+' = '+str(value)    
      #print "self."+name+' = '+str(value) 
  ##########################################################################################################    
  ###############################  convert results to new scale   ########################################## 
  ##########################################################################################################       
  def convertscale(self,pathandfile,pathfileout):
    if self.filetype=='SAS':
      target = open(pathandfile, 'r')

      #get index on which data is located for: var_scale_ref, var_dependent, and  var_fixed
      line1 = target.readline()
      line1 = line1.replace("\r\n","")
      line1 = line1.replace("\n","")      
      header = str.split(line1,',')
      
      var_scale_ref_index =  header.index(self.var_scale_ref)
      
      var_dependent_index = []
      for var in self.var_dependent:
        var_dependent_index.append(header.index(var))
      
      var_fixed_index = []  
      for var in self.var_fixed:
        var_fixed_index.append(header.index(var))

      #creat array to save data
      VSREFarray = []
      count=0
      for var in self.var_dependent:
        exec "VDEParray"+str(count)+" = []"
        count+=1
        
      count=0  
      for var in self.var_fixed:
        exec "VFIXarray"+str(count)+" = []"        
        count+=1
      '''
      internal managing:
        1 input data to lists
        2 reshape data for each Vd
        3 generate data using given Vg scale
        4 reshape to the original format
        5 write output file
      '''        
      #1-input data to arrays  
      flag2ndline = 1
      for line in target:
        #print line
        line = line.replace("\r\n","")
        line = line.replace("\n","")          
        linesplit = str.split(line,',')

        
        if flag2ndline==1:
          deviceinfo = line
          flag2ndline = 0
        VSREFarray.append(float(linesplit[var_scale_ref_index+1]))
        
        count=0
        for index in var_dependent_index:
          exec "VDEParray"+str(count)+".append(float(linesplit[index+1]))"
          count+=1
          
        count=0  
        for index in var_fixed_index:
          exec "VFIXarray"+str(count)+".append(float(linesplit[index+1]))"        
          count+=1
          
      target.close()
      #transform to numpy array    
      VSREFarray = np.array(VSREFarray)
        
      count=0
      for index in var_dependent_index:
        exec "VDEParray"+str(count)+"= np.array(VDEParray"+str(count)+")"
        count+=1
          
      count=0  
      for index in var_fixed_index:
        exec "VFIXarray"+str(count)+"= np.array(VFIXarray"+str(count)+")"       
        count+=1          
      
      #2 - reshape arrays  
      lengthalldata = len(VSREFarray)
      lengthsingleset = len(np.unique(VSREFarray)) 
      numberofsets = int(lengthalldata/lengthsingleset)

      #3 generate data using given Vg scale, 4 reshape to the original format          
      #3.a) generate reference array
      VSREFarraynewref = np.linspace(self.var_scale_ref_init,self.var_scale_ref_final,round((self.var_scale_ref_final-self.var_scale_ref_init)/self.var_scale_ref_delta))
      #VSREFarraynew = VSREFarraynew.reshape(len(VSREFarraynew),1)
      VSREFarraynew =  np.tile(VSREFarraynewref, numberofsets) 
      
      #3.b) generate array for dependent variables
      count=0
      VDEParraynew = np.array([])
      #this for is for each variable
      for index in var_dependent_index:
        countcolum=0
        #this for is for each set of given variable
        while countcolum<(numberofsets):
          stringexec =  "y = VDEParray"+str(count)+"[countcolum*lengthsingleset:(countcolum+1)*lengthsingleset-1]"
          #print stringexec
          exec stringexec
          tck = interpolate.splrep(VSREFarray[0:lengthsingleset-1], y, s=0)
          VDEParraynew = np.concatenate((VDEParraynew,interpolate.splev(VSREFarraynewref, tck, der=0)))
          countcolum+=1
        count+=1
      VDEParraynew = np.array(VDEParraynew)

      #3.c) generate fix array
      VFIXarraynew = np.array([])
      count=0  
      for index in var_fixed_index:
        exec "y = VFIXarray"+str(count)
        y = np.unique(y)
        VFIXarraynew = np.concatenate((VFIXarraynew,np.repeat(y.reshape(1,len(y)),len(VSREFarraynewref))))
        print y
        count+=1 
      
      ALLNEWarray = np.concatenate((VSREFarraynew , VDEParraynew , VFIXarraynew))
      ALLNEWarray = ALLNEWarray.reshape(len([self.var_scale_ref])+len(self.var_dependent)+len(self.var_fixed),len(ALLNEWarray)/(len([self.var_scale_ref])+len(self.var_dependent)+len(self.var_fixed)))

      #5 write output file
      indexinnew  = range(1+len(var_dependent_index)+len(var_fixed_index))
      indextoprint = [var_scale_ref_index]+var_dependent_index + var_fixed_index
      #namestoprint = [self.var_scale_ref]+self.var_dependent+self.var_fixed
      deviceinfolist = str.split(deviceinfo,',')      
      deviceinfolist = ','.join(deviceinfolist[0:-len(indextoprint)])
      indextoprint = [x for (y,x) in sorted(zip(indextoprint,indexinnew))]
      
      fileout = open(pathfileout, 'w')
      fileout.write(line1)  
            
      count=0
      while (count<(len(VSREFarraynew))):
        stringtoprint = deviceinfolist
        for index in indextoprint:
          stringtoprint +=','+str(ALLNEWarray[index,count]) 
        stringtoprint +='\n'
        fileout.write(stringtoprint) 
        count+=1
      fileout.close()
      #line1
      
  ##########################################################################################################    
  ##################################            plot              ########################################## 
  ########################################################################################################## 
  def plotfiledata(self,pathandfile,xstring,ystring,fignumber):
    #this function open a file pathandfile and plot the columns with xstring and ystring string header
    flagprint = 0
    if self.filetype=='':
      target = open( pathandfile, 'r')
      header = str.split(target.readline())
      if len(header)>0:
        flagprint = 1
        xindex =  header.index(xstring)
        yindex =  header.index(ystring)
      
        datalist = loadtxt(pathandfile,skiprows = 1)

        xarray = datalist[:,xindex]
        yarray = datalist[:,yindex]
      
      #this is to identify index how to re-shape matrix for right plotting
      numberelement = 0
      numberelementaux = len(np.unique(xarray)) 
      numberelementmaxpossible = len(xarray)
      if( (np.mod(len(xarray),numberelementaux)==0) and ((numberelementmaxpossible-numberelementaux)>0) ):
        numberelement = numberelementaux;
        elementpercylce = findelementpercylce(xarray)*numberelement
      
      if (numberelement==0):
        numberelement = numberelementmaxpossible
        elementpercylce = numberelement
      
      #reshape matrix to plot lines
      xarray = rearrangearray(xarray,elementpercylce,numberelement)
      yarray = rearrangearray(yarray,elementpercylce,numberelement)        
      target.close()
      
    #SAS format, for intership  
    if self.filetype=='SAS':
      flagprint = 1
      target = open(pathandfile, 'r')

      #datalist = loadtxt(pathdatafile,skiprows = 1)
      header = str.split(target.readline(),',')
      #print header
      xindex =  header.index(xstring)
      yindex =  header.index(ystring)
      #print xindex,yindex
      xarray = []
      yarray = []

      for line in target:
        #print line
        linesplit = str.split(line,',')
        #print linesplit
        if (linesplit[xindex+1]!='noData') and (linesplit[yindex+1]!='noData'):
          xarray.append(float(linesplit[xindex+1]))
          yarray.append(float(linesplit[yindex+1]))     
      target.close()
      
    if flagprint==1:  

      #plot
      plt.figure(fignumber)

      #log scale check  
      #yarray = abs(yarray)
      if self.ylogflag==1:
        yarray = abs(yarray) 
              
      #plot variable or its derivatives: TODO: it plot derivate with respect to x-axis, update derivative with respect to any variable
      if (self.derivativeorder<1):  
        if self.color=='':    
          plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize  )
        else:
          plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize, color=self.color  )
      else :
        K = K_generator(xarray[:,0],self.derivativeorder) 
        if self.color=='':    
          plt.plot( xarray, K*yarray, self.symbol, lw=self.lw,markersize=self.markersize)
        else:
          plt.plot( xarray, K*yarray, self.symbol, lw=self.lw,markersize=self.markersize, color=self.color  )
      
      #log scale check  
      if self.ylogflag==1:
        ax = plt.gca()
        ax.set_yscale('log')
      if self.xlogflag==1:
        ax = plt.gca()
        ax.set_xscale('log')        
    
    #x and y axis label          
    
      ax = plt.gca()
      ax.set_xlabel(xstring)
      if (self.derivativeorder<1):
        ax.set_ylabel(ystring)  
      else:
        ax.set_ylabel('d^'+str(self.derivativeorder)+' '+ystring+'/d'+xstring+'^'+str(self.derivativeorder))
        
  ##########################################################################################################    
  ###############################  add data to measurement analysis######################################### 
  ##########################################################################################################  
  #technology,wafer,site,macro,device,test,temperature,Vd`,Vg`,Ig,Id,Is
  #RsearchAtypical,KL35LC12-S01,18,16,35LC_GL1_M04,FET05,Idgs@Vg51d3,25, -5.000e-02,-5.000e-01,-7.9736e-11, 1.0170e-12, 4.9877e-11  
  
  '''
  1) identify input parameter names and order in head: ex, technology,wafer,site,macro,device,test,temperature,Vd`,Vg`,Ig,Id,Is 
  
  2) check if device exists: technology,wafer,site,macro,device
  
  3) chec if test+temperatu exists
  
  4) add data, until new device is found, repeat
  
  '''      
  def adddata(self,pathandfile):
    
    if self.filetype=='SAS':
      
      target = open(pathandfile, 'r')  
      state = -1 #reading header: technology,wafer,site,macro,device,test,temperature,Vd`,Vg`,Ig,Id,Is 
      countline=-1
      currentdevice = 'non'
      for line in target:
        countline+=1
        #print line.index('technology')
        #########################################################for reading head line
        if ('technology' in line):
          #print "A header has been found"
          line = line.replace("\r\n","")
          line = line.replace("\n","")            
          header = str.split(line,',')

          #obtain index for device tags: technology,wafer,site,macro,device,test,temperature
          device_tags_index = []
          for tag in self.device_tags:  
            device_tags_index.append(header.index(tag))
            
          #obtain index for the rest of variables
          measurement_parameters_index = []
          measurement_parameters_name = []
          count = 0
          while count<len(header):
            if not count in device_tags_index:
              measurement_parameters_index.append(count)
              measurement_parameters_name.append(header[count])
            count+=1
          #print measurement_parameters_index 
          #print measurement_parameters_name         
          state = 0
          currentdevice = ''
        ####################################################for reading rest of lines, and check if new device or test  
        if not 'technology' in line: 
          line = line.replace("\r\n","")
          line = line.replace("\n","")          
          header = str.split(line,',')
         
          if (currentdevice != ",".join(header[0:(self.device_tags.index('temperature')+2)])):
            print "\nNew measurement found"
            #print "currentdevice" + currentdevice
            print "New Device Data: " + ",".join(header[0:(self.device_tags.index('temperature')+2)]) 
            state = 1
        ######################################################input data to arrays in dictionaries  
        if state==1:
          state = 2
          line = line.replace("\r\n","")
          line = line.replace("\n","")          
          header = str.split(line,',')
          if (currentdevice != ",".join(header[0:(self.device_tags.index('temperature')+2)])):
            currentdevice = ",".join(header[0:(self.device_tags.index('temperature')+2)] )       
            namedevice= header[self.device_tags.index('technology')+1]
            namedevice+= ','+header[self.device_tags.index('wafer')+1]
            namedevice+= ','+header[self.device_tags.index('site')+1]
            namedevice+= ','+header[self.device_tags.index('macro')+1]
            namedevice+= ','+header[self.device_tags.index('device')+1]
            
            #check if device exist (check for key in dictionary
            #print namedevice
            if not namedevice in self.devices:
              print "Creating device named "  +namedevice
              self.devices[namedevice]  = {'technology':header[self.device_tags.index('technology')+1],'wafer':header[self.device_tags.index('wafer')+1],'site':header[self.device_tags.index('site')+1],'macro':header[self.device_tags.index('macro')+1],'device':header[self.device_tags.index('device')+1],'temperatures':[],'tests':[]} 
            
            #check if temperature exists 
            temperature = header[self.device_tags.index('temperature')+1]
            if not temperature in self.devices[namedevice]:
              print "Device named "+namedevice+ ", creates a new temperature data: "+ temperature  + " degrees"
              self.devices[namedevice][temperature]  = {}           
              self.devices[namedevice]['temperatures'].append(temperature)
            #check if test exists
            testname = header[self.device_tags.index('test')+1]
            if not testname in self.devices[namedevice][temperature]:
              print "Device named "+namedevice+ ", with temperature data: "+ temperature  + " degrees, creates new test: " + testname
              self.devices[namedevice][temperature][testname]  = self.device_template
              self.devices[namedevice]['tests'].append(testname)
              for varname in measurement_parameters_name:
                self.devices[namedevice][temperature][testname][varname] = [] 
              
        #####################################################add data to arrays in dictionaries
        if state == 2:
          line = line.replace("\r\n","")
          line = line.replace("\n","")
          linedata = str.split(line,',')

          
          count=0
          for varname in measurement_parameters_name:
            value = float(linedata[measurement_parameters_index[count]+1])
            self.devices[namedevice][temperature][testname][varname].append(value) 
            count+=1
      #print self.devices[namedevice][temperature][testname]['Vg`']      
      target.close()  
  ##########################################################################################################    
  ###############################  plot devices ######################################### 
  ##########################################################################################################        
  def plotdevices(self,fignumber):
    #plot method for data that is already inputed
    for technology in self.plot_technology:
      for wafer in self.plot_wafer:
        for site in self.plot_site:
          for macro in self.plot_macro:
            for device in self.plot_device:
              for temperature in self.plot_temperature:
                for test in self.plot_test:
                  namedevice = technology+','+wafer+','+site+','+macro+','+device
                  if namedevice in self.devices:
                    if temperature in self.devices[namedevice]:
                      if test in self.devices[namedevice][temperature]:
                        xarray = self.devices[namedevice][temperature][test][self.plot_x_variable]
                        #print xarray
                        for ystring in self.plot_y_variables:
                          #print "ploting: " + ystring
                          yarray = self.devices[namedevice][temperature][test][ystring]  
                          plt.figure(fignumber)

                          #log scale check  
                          if self.ylogflag==1:
                            yarray = abs(yarray) 
                                  
                          #plot variable or its derivatives: TODO: it plot derivate with respect to x-axis, update derivative with respect to any variable
                          if (self.derivativeorder<1):  
                            if self.color=='':    
                              plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize  )
                            else:
                              plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize, color=self.color  )
                          else :
                            K = K_generator(xarray[:,0],self.derivativeorder) 
                            if self.color=='':    
                              plt.plot( xarray, K*yarray, self.symbol, lw=self.lw,markersize=self.markersize)
                            else:
                              plt.plot( xarray, K*yarray, self.symbol, lw=self.lw,markersize=self.markersize, color=self.color  )
      #log scale check  
    if self.ylogflag==1:
      ax = plt.gca()
      ax.set_yscale('log')
    if self.xlogflag==1:
      ax = plt.gca()
      ax.set_xscale('log')        
    
    #x and y axis label          
    
    ax = plt.gca()
    ax.set_xlabel(self.plot_x_variable)
    ylegend = ''
    for ystring in self.plot_y_variables:
      ylegend+=ystring+' '
    if (self.derivativeorder<1):
      ax.set_ylabel(ylegend)  
    else:
      ax.set_ylabel('d^'+str(self.derivativeorder)+' '+ylegend+'/d'+self.plot_x_variable+'^'+str(self.derivativeorder))       
  ##########################################################################################################    
  ###############################  add all data from folder ######################################### 
  ####################################################################################################                   
  def addalldatainfolder(self,folder):
    data_files = [(x[0], x[2]) for x in os.walk(folder)]
    for datafile in data_files[0][1]:
      print "\n Adding: "+datafile
      self.adddata(folder+datafile)
                
