#this class is to handle data obtained from measurements
#Juan Pablo Duarte
#########################################################################
import os
import matplotlib.pyplot as plt
from numpy import loadtxt
import numpy as np
from scipy import interpolate

import supportfunctions as sp
import devicecharacterization as devc

#########################################################################
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
    self.device_info = {}
    self.plot_all_together = 1
    self.plot_characterization_fit=0
    self.plot_characterization_save = 0
    self.plot_characterization_file_out = ''
    self.delta_Lg = 0
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
        elementpercylce = sp.findelementpercylce(xarray)*numberelement
      
      if (numberelement==0):
        numberelement = numberelementmaxpossible
        elementpercylce = numberelement
      
      #reshape matrix to plot lines
      xarray = sp.rearrangearray(xarray,elementpercylce,numberelement)
      yarray = sp.rearrangearray(yarray,elementpercylce,numberelement)        
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
        K = sp.K_generator(xarray[:,0],self.derivativeorder) 
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
      countdata=0
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
        if (not 'technology' in line) and (state!=-1): 
          line = line.replace("\r\n","")
          line = line.replace("\n","")          
          header = str.split(line,',')
         
          if (currentdevice != ",".join(header[0:(self.device_tags.index('temperature')+2)])):
            #if countdata>1:
              #print self.devices[namedevice][temperature]['Idgsx@Vg126d2']['Vg`']
            print "Previous data size: "+str(countdata)
            countdata=0
            print "\nNew measurement found"
            print "New device" + ",".join(header[0:(self.device_tags.index('temperature')+2)])
            print "New Device Data: " + ",".join(header[0:(self.device_tags.index('temperature')+2)]) 
            
            state = 1
        ######################################################input data to arrays in dictionaries  
        if state==1:
          countdata+=1
          state = 2
          line = line.replace("\r\n","")
          line = line.replace("\n","")          
          header = str.split(line,',')
          print header
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
              
              self.devices[namedevice]['Wdes'] = self.device_info[namedevice]['Wdes']
              self.devices[namedevice]['Ldes'] = self.device_info[namedevice]['Ldes']
              self.devices[namedevice]['Leff'] = self.device_info[namedevice]['Leff'] 
            
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
              self.devices[namedevice][temperature][testname]  = {}
              self.devices[namedevice]['tests'].append(testname)
              for varname in measurement_parameters_name:
                self.devices[namedevice][temperature][testname][varname] = [] 
              
        #####################################################add data to arrays in dictionaries
        if state == 2:
          countdata+=1
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
  ###############################  add data info of wafer die######################################### 
  ##########################################################################################################        
  def add_die_info(self,pathandfile):
    
    if self.filetype=='SAS':
      
      target = open(pathandfile, 'r')  
      state = -1 #reading header: technology,wafer,site,macro,device,test,temperature,Vd`,Vg`,Ig,Id,Is 

      for line in target:
        #print line.index('technology')
        #########################################################for reading head line
        if state==-1:
          #print "A header has been found"
          line = line.replace("\t"," ")
          line = line.replace("\r\n","")
          line = line.replace("\n","")            
          headernames = str.split(line)
          #print "Header: "+line     
          state = 0
        else:
          line = line.replace("\t"," ")
          line = line.replace("\r\n","")
          line = line.replace("\n","")          
          header = str.split(line)
          if True:
            #DevNum Wafer SiteX SiteY Macro Device  Wdes  Ldes
            #print headernames
            namedevice= header[headernames.index('Wafer')]
            namedevice+= ','+header[headernames.index('SiteX')]
            namedevice+= ','+header[headernames.index('SiteY')]
            namedevice+= ','+header[headernames.index('Macro')]
            namedevice+= ','+header[headernames.index('Device')]
            
            #check if device exist (check for key in dictionary
            #print namedevice
            if not namedevice in self.device_info:
              #print "Creating device named "  +namedevice
              self.device_info[namedevice]  = {'DevNum':header[headernames.index('DevNum')],'Wafer':header[headernames.index('Wafer')],'SiteX':header[headernames.index('SiteX')],'SiteY':header[headernames.index('SiteY')],'Macro':header[headernames.index('Macro')],'Device':header[headernames.index('Device')],'Wdes':float(header[headernames.index('Wdes')]),'Ldes':float(header[headernames.index('Ldes')]),'Leff':float(header[headernames.index('Ldes')])-self.delta_Lg} #TODO: do not hard code 0.04 for Leff calculation
    
      target.close()        
  ##########################################################################################################    
  ###############################  plot devices ######################################### 
  ##########################################################################################################        
  def plotdevices(self,fignumber):
    #plot method for data that is already inputed

    if self.plot_technology == 'all':
      technology_array = self.allinforeturn('technology')
    else:
      technology_array = self.plot_technology
    for technology in technology_array:    
      
      if (self.plot_wafer == 'all'):
        wafer_array = self.allinforeturn('wafer')
      else:
        wafer_array = self.plot_wafer    
      for wafer in wafer_array:
      
        if self.plot_site == 'all':
          site_array = self.allinforeturn('site')
        else:
          site_array = self.plot_site       
        for site in site_array:
        
          if self.plot_macro == 'all':
            macro_array = self.allinforeturn('macro')
          else:
            macro_array = self.plot_macro          
          for macro in macro_array:
          
            if self.plot_device == 'all':
              device_array = self.allinforeturn('device')
            else:
              device_array = self.plot_device             
            for device in device_array:
            
              if self.plot_temperature == 'all':
                temperature_array = self.allinforeturn('temperatures')
              else:
                temperature_array = self.plot_temperature             
              for temperature in temperature_array:
                for test in self.plot_test:
                  namedevice = technology+','+wafer+','+site+','+macro+','+device
                  if namedevice in self.devices:
                    if temperature in self.devices[namedevice]:
                      if test in self.devices[namedevice][temperature]:
                        xarray = self.devices[namedevice][temperature][test][self.plot_x_variable]
                        #print xarray
                        for ystring in self.plot_y_variables:
                          print "ploting: " + ystring
                          yarray = self.devices[namedevice][temperature][test][ystring]
                          if self.plot_W_normalization_flag==1:
                            yarray = np.array(yarray)/self.devices[namedevice]['Wdes']
                          #########################update array for given limiting conditions
                          legendlimit = []
                          count=0
                          for limitinfo in self.plot_parameter_limits:
                            legendlimit.append([])
                            count+=1

                          yaux = []
                          xaux = []
                          count=0
                          flagall = 0
                          for yvalue in yarray:
                            flag=1
                            countlimit=0
                            for limitinfo in self.plot_parameter_limits:
                              if limitinfo[0] in self.devices[namedevice][temperature][test]:
                                valueaux = self.devices[namedevice][temperature][test][limitinfo[0]][count]
                              elif limitinfo[0] in self.devices[namedevice]:
                                valueaux = self.devices[namedevice][limitinfo[0]]
                              else:
                                print "ERROR in limit"
                                
                              if len(limitinfo[1])==1:
                                maxvalue = limitinfo[1][0]+1e-14
                                minvalue = limitinfo[1][0]-1e-14
                              else:
                                maxvalue = max(limitinfo[1])
                                minvalue = min(limitinfo[1])                             
                              if not ((valueaux<maxvalue) and (valueaux>minvalue)):
                                flag=0 
                              countlimit+=1
                            if flag==1:
                              flagall=1
                              yaux.append(yvalue)
                              xaux.append(xarray[count])
                              countaddlegen=0
                              for limitinfo in self.plot_parameter_limits:
                                if limitinfo[0] in self.devices[namedevice][temperature][test]:
                                  valueaux = self.devices[namedevice][temperature][test][limitinfo[0]][count]
                                elif limitinfo[0] in self.devices[namedevice]:
                                  valueaux = self.devices[namedevice][limitinfo[0]]
                                else:
                                  print "ERROR in limit"                                
                                legendlimit[countaddlegen].append(valueaux)
                                countaddlegen+=1
                              
                            count+=1  
                          
                          if flagall==1:
                            yarray =  yaux
                            xarray =  xaux   
                            
                            #############################################check if one plot per data or all together
                            if self.plot_all_together ==1: 
                              plt.figure(fignumber)
                            else:
                              plt.figure(fignumber)
                              fignumber+=1

                            #log scale check  
                            if self.ylogflag==1:
                              yarray = abs(np.array(yarray)) 
                            
                            ###################################create legend
                            namelegend = ''
                            for stringaux in self.plot_legend_names:   
                              if (stringaux == 'bias'):
                                biasstring = ''
                                countleg=0
                                for limitinfo in self.plot_parameter_limits:
                                  biasstring = biasstring + limitinfo[0] + '='+','.join(str(x) for x in list(set(legendlimit[countleg])) )+', '
                                  countleg+=1
                                namelegend=namelegend+biasstring
                              elif stringaux in ['Leff','Wdes','Ldes']:
                                namelegend=namelegend+stringaux +':'+str(self.devices[namedevice][stringaux]) + ', '
                              else:
                                namelegend=namelegend+locals()[stringaux] + ', '
                            namelegend = namelegend[0:-2]
                            ##################################create title  
                            nametitle = ''
                            for stringaux in self.plot_title_names : 
                              nametitle=nametitle+locals()[stringaux] + ','
                            ax = plt.gca()
                            ax.set_title(nametitle)                                                         
                            
                            #plot variable or its derivatives: TODO: it plot derivate with respect to x-axis, update derivative with respect to any variable
                            if (self.derivativeorder<1):  
                              if self.color=='':    
                                plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize  , label=namelegend)
                              else:
                                plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize, color=self.color  , label=namelegend)
                            else :
                              K = sp.K_generator(xarray[:,0],self.derivativeorder) 
                              if self.color=='':    
                                plt.plot( xarray, K*yarray, self.symbol, lw=self.lw,markersize=self.markersize, label=namelegend)
                              else:
                                plt.plot( xarray, K*yarray, self.symbol, lw=self.lw,markersize=self.markersize, color=self.color  , label=namelegend)
      #log scale check  
    if self.ylogflag==1:
      ax = plt.gca()
      ax.set_yscale('log')
      plt.legend(loc='lower right')#, bbox_to_anchor=(1, 0.5))
    else:
      plt.legend(loc='upper left')#, bbox_to_anchor=(1, 0.5))
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
      ylegend = ylegend 
    else:
      ylegend = 'd^'+str(self.derivativeorder)+' '+ylegend+'/d'+self.plot_x_variable+'^'+str(self.derivativeorder)
    if self.plot_W_normalization_flag==1:
      ylegend=ylegend+'/um'
        
    ax.set_ylabel(ylegend)   
  ##########################################################################################################    
  ###############################  add all data from folder ######################################### 
  ####################################################################################################                   
  def addalldatainfolder(self,folder):
    data_files = [(x[0], x[2]) for x in os.walk(folder)]
    count=1
    for datafile in data_files[0][1]:
      print "\n Adding folder number " + str(count)+': '+datafile
      self.adddata(folder+datafile)
      count+=1
      
  ##########################################################################################################    
  ###############################  calculate device electrical quantities ################################## 
  ##########################################################################################################         
  def device_characterization(self):
    #Vth for different Vds, Tempt
    print "...running electrical characterization of devices..."
    '''
    ['Ix', 'Is', 'Vg`', 'Vd`', 'Id', 'Ig']
M1.updateparameter('vth_testname' ,'Idgsx@Vg126d2')
M1.updateparameter('vth_biasreference' , 'Vg`')
M1.updateparameter('vth_method' , 'constant_current')
M1.updateparameter('vth_biasfixed' , ['Vd`'])
M1.updateparameter('vth_current_level' , 300e-9)
    '''
    
    #Vth, and DIBL calculation
    for device in self.devices.keys():
      for temperature in self.devices[device]['temperatures']:
        Vref = np.array(self.devices[device][temperature][self.vth_testname][self.vth_biasreference])
        Vfixed = np.array(self.devices[device][temperature][self.vth_testname][self.vth_biasfixed[0]])#TODO: include more biases, for loop
        Ids = np.array(self.devices[device][temperature][self.vth_testname]['Id'])
        
        numbervths = len(np.unique(Vfixed))
        lengthalldata = len(Vref)
        Vdsunique = np.unique(Vfixed)
        
        Vth = []
        SS  = []
        Gmmax = []
        count=0
        Ionmax = []
        Ronmax = []
        while (count<numbervths):
          idsref = self.vth_current_level*self.devices[device]['Wdes']/self.devices[device]['Leff']
          Vrefaux = Vref[count*lengthalldata/numbervths:(count+1)*lengthalldata/numbervths]
          Idsaux = Ids[count*lengthalldata/numbervths:(count+1)*lengthalldata/numbervths]
          #print "Vth calculations"
          vthaux = devc.findvth(Vrefaux,Idsaux,idsref)
          Vth.append(vthaux)
          SS.append(devc.findss(Vrefaux,Idsaux,vthaux,self.ss_vthwindow,self.ss_fitwindow))
          Gmmax.append(devc.findGmax(Vrefaux,Idsaux)/self.devices[device]['Wdes'])
          Ionmax.append(max(Idsaux)/self.devices[device]['Wdes'])
          Ronmax.append(Vdsunique[count]/(max(Idsaux)/self.devices[device]['Wdes']) )
          
          count+=1
        DIBL = (max(Vth)-min(Vth))/(max(Vfixed)-min(Vfixed))
        self.devices[device][temperature]['characterization'] = {'vth':Vth,'bias_fixed': np.unique(Vfixed), 'dibl': DIBL, 'ss':SS, 'gmmax': Gmmax, 'ionmax':Ionmax,'ronmax':Ronmax}
    '''
    M1.updateparameter('gmmax_method' , 'ioffref')
M1.updateparameter('ion_method' , 'ioffref')
M1.updateparameter('ron_method' , 'contant_overdrive')
M1.updateparameter('ron_overdrive' , 0.7)
M1.updateparameter('vdd_ref' , 0.5)
M1.updateparameter('ioff_ref' , 0.5)
    '''
    
    for device in self.devices.keys():
      for temperature in self.devices[device]['temperatures']:
        Vref = np.array(self.devices[device][temperature][self.vth_testname][self.vth_biasreference])
        Vfixed = np.array(self.devices[device][temperature][self.vth_testname][self.vth_biasfixed[0]])#TODO: include more biases, for loop
        Ids = np.array(self.devices[device][temperature][self.vth_testname]['Id'])
        
        numbervths = len(np.unique(Vfixed))
        lengthalldata = len(Vref)
        Vdsunique = np.unique(Vfixed)
        
        Vthlin = max(self.devices[device][temperature]['characterization']['vth'])
        Vgon = Vthlin+self.ron_overdrive
        VgIoff = Vthlin+self.vgs_off
        VgIon  = VgIoff+self.vgs_off
        
        Ron = []
        Ion = []
        Ioff = []
        count=0
        while (count<numbervths):
          idsref = self.vth_current_level*self.devices[device]['Wdes']/self.devices[device]['Leff']
          Vrefaux = Vref[count*lengthalldata/numbervths:(count+1)*lengthalldata/numbervths]
          Idsaux = Ids[count*lengthalldata/numbervths:(count+1)*lengthalldata/numbervths]
          Ron.append(devc.findRon(Vrefaux,Idsaux,Vgon,Vdsunique[count])*self.devices[device]['Wdes'])
          SS = self.devices[device][temperature]['characterization'] ['ss'][count]
          Ion.append(devc.findIds(Vrefaux,Idsaux,VgIon,SS)/self.devices[device]['Wdes'])
          Ioff.append(devc.findIds(Vrefaux,Idsaux,VgIoff,SS)/self.devices[device]['Wdes'])
          count+=1
          
        self.devices[device][temperature]['characterization']['ron'] = Ron
        self.devices[device][temperature]['characterization']['ion'] = Ion
        self.devices[device][temperature]['characterization']['ioff'] = Ioff
        
        #self.devices[device][temperature]['characterization']['SS'] = SS
  ##########################################################################################################    
  ###############################  extract all the tech,macro,site,etc, from the dictionary################# 
  ##########################################################################################################    
  def allinforeturn(self,tag):
  
    tagtoreturn = []
    for device in self.devices.keys():
      tagtoreturn.append(self.devices[device][tag])            
    return list(set(tagtoreturn))
       
  def plotcharacterization(self,fignumber,typeplot,xvariable,extralegend):
  
    if True:
    #Vth lin vs Lg
      xarray = []
      yarray = []
      for device in self.devices.keys():
        for temperature in self.devices[device]['temperatures']:
          #self.devices[device][temperature]['characterization'] = {'Vth':Vth,'bias_fixed': np.unique(Vfixed), 'DIBL': DIBL, 'SS':SS}
          indexvth = sp.find_closest(self.devices[device][temperature]['characterization']['bias_fixed'], self.plot_characterization_vdref)
          if typeplot.lower()=='dibl':
            yarray.append(self.devices[device][temperature]['characterization'][typeplot.lower()])
          else:
            yarray.append(self.devices[device][temperature]['characterization'][typeplot.lower()][indexvth])
          if xvariable in self.devices[device]:
            xarray.append(self.devices[device][xvariable])
          else:
            xarray.append(self.devices[device][temperature]['characterization'][xvariable.lower()][indexvth])
      plt.figure(fignumber) 
      if typeplot.lower()=='dibl':
        namelegend = typeplot
      else:
        namelegend = typeplot+'@Vds:'+str(self.plot_characterization_vdref) 
      namelegend=namelegend+extralegend
      if self.ylogflag==1:
        yarray = abs(np.array(yarray))
      if self.xlogflag==1:
        xarray = abs(np.array(xarray))  
              
      if self.color=='':    
        plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize  , label=namelegend)
      else:
        plt.plot( xarray, yarray, self.symbol, lw=self.lw,markersize=self.markersize, color=self.color  , label=namelegend)
      #log scale check  
      if (self.plot_characterization_fit >0):
        z = np.polyfit(np.array(xarray), np.array(yarray), self.plot_characterization_fit)
        print "Poly Fit: "
        print z
        p = np.poly1d(z)
        xaux = np.linspace(min(xarray),max(xarray),100)
        if self.color=='': 
          plt.plot( xaux, p(xaux), '-', lw=self.lw,markersize=self.markersize  , label=namelegend+'fit')
        else:
          plt.plot( xaux, p(xaux), '-', lw=self.lw,markersize=self.markersize, color=self.color  , label=namelegend+'fit')
      
    if self.ylogflag==1:
      ax = plt.gca()
      ax.set_yscale('log')
      plt.legend(loc='lower right')#, bbox_to_anchor=(1, 0.5))
    else:
      plt.legend(loc='upper left')#, bbox_to_anchor=(1, 0.5))
    if self.xlogflag==1:
      ax = plt.gca()
      ax.set_xscale('log')        
    
    #x and y axis label          
    
    ax = plt.gca()
    ax.set_xlabel(xvariable)      
    ax.set_ylabel(typeplot) 
    
    if (self.plot_characterization_save==1):
      fileresult = open(self.plot_characterization_file_out, 'w') 
      fileresult.write(xvariable+' '+typeplot+' \n')
      #this save to text file by evaluating model
      i=0
      column=len(xarray)
      rown=len(xarray)
      while (i<rown):
        stringtowrite = str(xarray[i])+' '+ str(yarray[i])
        fileresult.write(stringtowrite+'\n') 
        i+=1 
      fileresult.close()
