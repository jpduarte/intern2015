#script example to input data from several devices, and plot it
#Juan Pablo Duarte

rootfolder = '/home/juan/research/intern2015'

#indicate path for folders containing required classes
import sys
sys.path.insert(0, rootfolder+'/classesintern')

#import class for python simulation
import matplotlib.pyplot as plt
import measurementtests as mt
import supportfunctions as sp


figurename = '50kvs300k_'
###################################definition of class mt
M1 = mt.mt()

###################################input data for different devices
pathfolder = '/home/juan/research/intern/data/'

M1.updateparameter('device_tags',['technology','wafer','site','macro','device','test','temperature'])

M1.add_die_info('/home/juan/research/intern/SC_III-V/SA35FET19.txt')
#add entire
M1.addalldatainfolder('/home/juan/research/intern/SC_III-V/50K/')

#######################################plot set up
M1.updateparameter('plot_technology','all')
M1.updateparameter('plot_wafer','all')
M1.updateparameter('plot_site','all')
M1.updateparameter('plot_macro','all')#width,35fet_GF2_FetPS_10
M1.updateparameter('plot_device','all')#->gate length
M1.updateparameter('plot_test',['Idgs@Vg23d2'])
M1.updateparameter('plot_temperature',['25'])
M1.updateparameter('plot_y_variables',['Id'])
M1.updateparameter('plot_x_variable','Vg`')
M1.updateparameter('plot_legend_names',['Leff','Wdes','bias'])#
M1.updateparameter('plot_title_names',['technology','wafer','site','macro'])
M1.updateparameter('plot_parameter_limits',[ ['Vd`',[0.05]] ])#, ['Leff',[0.2,0.3]]
M1.updateparameter('plot_W_normalization_flag',1)#
######################################plot run
M1.updateparameter('plot_all_together' , 1)
#M1.plotdevices(1)

M1.updateparameter('ylogflag' , 1)
#M1.plotdevices(2)

#######################################characterization set up
M1.updateparameter('vth_testname' ,'Idgsx@Vg126d2')
M1.updateparameter('vth_biasreference' , 'Vg`')
M1.updateparameter('vth_method' , 'constant_current')
M1.updateparameter('vth_biasfixed' , ['Vd`'])
M1.updateparameter('vth_current_level' , 300e-9)
M1.updateparameter('ss_method' , 'max_ss')
M1.updateparameter('ss_vthwindow' , 0.5)
M1.updateparameter('ss_fitwindow' , 0.2)
M1.updateparameter('gmmax_method' , 'ioffref')
M1.updateparameter('ion_method' , 'ioffref')
M1.updateparameter('ron_method' , 'constant_overdrive')
M1.updateparameter('ron_overdrive' , 0.7)
M1.updateparameter('vdd_ref' , 0.5)
M1.updateparameter('vgs_off' , -0.3)
#M1.updateparameter('ioff_ref' , 0.5)

M1.device_characterization()

M1.updateparameter('ylogflag' , 0)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'r')
M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.plotcharacterization(1,'Vth','Leff','-50K')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'r')
M1.plotcharacterization(2,'Vth','Leff','-50K')

M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'r')
M1.plotcharacterization(3,'DIBL','Leff','-50K')

M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('color' , 'r')
M1.updateparameter('symbol' , 's')
M1.plotcharacterization(4,'SS','Leff','-50K')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('symbol' , 'o')
M1.updateparameter('color' , 'r')
M1.plotcharacterization(5,'SS','Leff','-50K')
#0-2

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('xlogflag' , 0)
M1.updateparameter('color' , 'r')
M1.plotcharacterization(6,'Ron','Leff','-50K')
#0-1000

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'r')
M1.plotcharacterization(7,'Gmmax','Leff','-50K')

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'r')
M1.plotcharacterization(8,'Ion','Leff','-Juan')

########################################show plots
printflag=1
plothead = '300Kvs50K_'

###################################definition of class mt
M1 = mt.mt()

###################################input data for different devices
pathfolder = '/home/juan/research/intern/data/'

M1.updateparameter('device_tags',['technology','wafer','site','macro','device','test','temperature'])

M1.add_die_info('/home/juan/research/intern/SC_III-V/SA35FET19.txt')
#add entire
M1.addalldatainfolder('/home/juan/research/intern/SC_III-V/300Kv2/')

#######################################plot set up
M1.updateparameter('plot_technology','all')
M1.updateparameter('plot_wafer','all')
M1.updateparameter('plot_site','all')
M1.updateparameter('plot_macro','all')#width,35fet_GF2_FetPS_10
M1.updateparameter('plot_device','all')#->gate length
M1.updateparameter('plot_test',['Idgsx@Vg126d2'])
M1.updateparameter('plot_temperature',['25'])
M1.updateparameter('plot_y_variables',['Id'])
M1.updateparameter('plot_x_variable','Vg`')
M1.updateparameter('plot_legend_names',['Leff','Wdes','bias'])#
M1.updateparameter('plot_title_names',['technology','wafer','site','macro'])
M1.updateparameter('plot_parameter_limits',[ ['Vd`',[0.05]] ])#, ['Leff',[0.2,0.3]]
M1.updateparameter('plot_W_normalization_flag',1)#
######################################plot run
M1.updateparameter('plot_all_together' , 1)
#M1.plotdevices(1)

M1.updateparameter('ylogflag' , 1)
#M1.plotdevices(2)

#######################################characterization set up
M1.updateparameter('vth_testname' ,'Idgsx@Vg126d2')
M1.updateparameter('vth_biasreference' , 'Vg`')
M1.updateparameter('vth_method' , 'constant_current')
M1.updateparameter('vth_biasfixed' , ['Vd`'])
M1.updateparameter('vth_current_level' , 300e-9)
M1.updateparameter('ss_method' , 'max_ss')
M1.updateparameter('ss_vthwindow' , 0.5)
M1.updateparameter('ss_fitwindow' , 0.2)
M1.updateparameter('gmmax_method' , 'ioffref')
M1.updateparameter('ion_method' , 'ioffref')
M1.updateparameter('ron_method' , 'constant_overdrive')
M1.updateparameter('ron_overdrive' , 0.7)
M1.updateparameter('vdd_ref' , 0.5)
M1.updateparameter('vgs_off' , -0.3)
#M1.updateparameter('ioff_ref' , 0.5)

M1.device_characterization()

M1.updateparameter('ylogflag' , 0)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('color' , 'b')
M1.plotcharacterization(1,'Vth','Leff','-Juan')
plt.savefig(figurename+'vthlin.png', bbox_inches='tight')
if printflag==1:
  plt.savefig(plothead+'vthlin.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'b')
M1.plotcharacterization(2,'Vth','Leff','-Juan')
plt.savefig(figurename+'vthsat.png', bbox_inches='tight')
if printflag==1:
  plt.savefig(plothead+'vthsat.png', bbox_inches='tight')

M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'b')
M1.plotcharacterization(3,'DIBL','Leff','-Juan')
plt.savefig(figurename+'dibl.png', bbox_inches='tight')
if printflag==1:
  plt.savefig(plothead+'dibl.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('color' , 'b')
M1.updateparameter('symbol' , 's')
M1.plotcharacterization(4,'SS','Leff','-Juan')
axes = plt.gca()
axes.set_ylim([0,2])
plt.savefig(figurename+'sslin.png', bbox_inches='tight')
if printflag==1:
  plt.savefig(plothead+'sslin.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('color' , 'b')
M1.updateparameter('symbol' , 'o')
M1.plotcharacterization(5,'SS','Leff','-Juan')
axes = plt.gca()
axes.set_ylim([0,2])
plt.savefig(figurename+'sssat.png', bbox_inches='tight')
if printflag==1:
  plt.savefig(plothead+'sssat.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('xlogflag' , 0)
M1.updateparameter('color' , 'b')
M1.plotcharacterization(6,'Ron','Leff','-Juan')
axes = plt.gca()
axes.set_ylim([0,1000])
axes = plt.gca()
axes.set_xlim([0,0.5])
plt.savefig(figurename+'ron.png', bbox_inches='tight')
if printflag==1:
  plt.savefig(plothead+'ron.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'b')
M1.plotcharacterization(7,'Gmmax','Leff','-Juan')
axes = plt.gca()
axes.set_ylim([0,0.005])
plt.savefig(figurename+'gmax.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'b')
M1.plotcharacterization(8,'Ion','Leff','-Juan')
plt.savefig(figurename+'Ion.png', bbox_inches='tight')
if printflag==1:
  plt.savefig(plothead+'IonvsLeff.png', bbox_inches='tight')
########################################show plots

plt.show() 





