#script example to input data from several devices, and plot it
#Juan Pablo Duarte

rootfolder = '/home/juan/research/intern2015'

#indicate path for folders containing required classes
import sys
sys.path.insert(0, rootfolder+'/classesintern')

#import class for python simulation
import matplotlib.pyplot as plt
import measurementtests as mt

###################################definition of class mt
M1 = mt.mt()

###################################input data for different devices
pathfolder = '/home/juan/research/intern/data/'

M1.updateparameter('device_tags',['technology','wafer','site','macro','device','test','temperature'])

M1.add_die_info('/home/juan/research/intern/SC_III-V/SA35FET19.txt')
#add entire
M1.addalldatainfolder('/home/juan/research/intern/SC_III-V/300K/')

#######################################plot set up
M1.updateparameter('plot_technology',['KL35LC12-S01'])
M1.updateparameter('plot_wafer',['18'])
M1.updateparameter('plot_site',['16'])
M1.updateparameter('plot_macro',['35LC_GL1_M04'])
M1.updateparameter('plot_device',['FET07'])
M1.updateparameter('plot_test',['Idgs@Vg51d3'])
M1.updateparameter('plot_temperature',['25'])

M1.updateparameter('plot_y_variables',['Id'])
M1.updateparameter('plot_x_variable','Vg`')

count=1
for device in M1.devices.keys():
  print "Device count: " + str(count)
  print M1.devices[device]['Leff']
  print M1.devices[device]['temperatures']
  count+=1
######################################plot run
#M1.plotdevices(1)
#plt.show() 




