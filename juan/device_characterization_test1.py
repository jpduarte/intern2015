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
M1.updateparameter('plot_technology','all')
M1.updateparameter('plot_wafer','all')
M1.updateparameter('plot_site','all')
M1.updateparameter('plot_macro','all')#width,35fet_GF2_FetPS_10
M1.updateparameter('plot_device',['FET05'])#->gate length
M1.updateparameter('plot_test',['Idgsx@Vg126d2'])
M1.updateparameter('plot_temperature',['25'])
M1.updateparameter('plot_y_variables',['Id'])
M1.updateparameter('plot_x_variable','Vg`')
M1.updateparameter('plot_legend_names',['Leff','Wdes','bias'])#
M1.updateparameter('plot_title_names',['technology','wafer','site','macro'])
M1.updateparameter('plot_parameter_limits',[['Vd`',[0.5]],['Leff',[0.01,0.03]]])#
######################################plot run
M1.updateparameter('plot_all_together' , 1)
M1.plotdevices(1)

M1.updateparameter('ylogflag' , 1)
M1.plotdevices(2)

#######################################characterization set up
M1.updateparameter('vth_testname' ,'Idgsx@Vg126d2')
M1.updateparameter('vth_biasreference' , 'Vg`')
M1.updateparameter('vth_method' , 'constant_current')
M1.updateparameter('vth_biasfixed' , ['Vd`'])
M1.updateparameter('vth_current_level' , 300e-9)

M1.device_characterization()



########################################show plots
plt.show() 




