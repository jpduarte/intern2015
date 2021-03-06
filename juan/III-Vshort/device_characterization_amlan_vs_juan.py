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

###################################definition of class mt
M1 = mt.mt()

###################################input data for different devices
pathfolder = '/home/juan/research/intern/data/'

M1.updateparameter('device_tags',['technology','wafer','site','macro','device','test','temperature'])

M1.add_die_info('/home/juan/research/intern/SC_III-V/SA35FET19.txt')
#add entire
M1.addalldatainfolder('/home/juan/research/intern/SC_III-V/300Kamlan/')

#######################################plot set up
M1.updateparameter('plot_technology','all')
M1.updateparameter('plot_wafer','all')
M1.updateparameter('plot_site','all')
M1.updateparameter('plot_macro','all')#width,35fet_GF2_FetPS_10
M1.updateparameter('plot_device','all')#->gate length
M1.updateparameter('plot_test',['Idgs@Vg23d2'])
M1.updateparameter('plot_temperature',['300'])
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
M1.updateparameter('vth_testname' ,'Idgs@Vg23d2')
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
M1.plotcharacterization(1,'Vth','Leff','--Amlan')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'r')
M1.plotcharacterization(2,'Vth','Leff','-Amlan')

M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'r')
M1.plotcharacterization(3,'DIBL','Leff','-Amlan')

M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('color' , 'r')
M1.updateparameter('symbol' , 's')
M1.plotcharacterization(4,'SS','Leff','-Amlan')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('symbol' , 'o')
M1.updateparameter('color' , 'r')
M1.plotcharacterization(5,'SS','Leff','-Amlan')
#0-2

M1.updateparameter('plot_characterization_fit' , 1)
M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('xlogflag' , 0)
M1.updateparameter('color' , 'r')
M1.plotcharacterization(6,'Ron','Leff','-Amlan')
#0-1000

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'r')
M1.plotcharacterization(7,'Gmmax','Leff','-Amlan')

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('ylogflag' , 1)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'r')
M1.plotcharacterization(8,'Ioff','Ion','-Amlan')

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('ylogflag' , 1)
M1.updateparameter('xlogflag' , 0)
M1.updateparameter('color' , 'r')
M1.plotcharacterization(9,'Iofffloor','Leff','-Amlan')
########################################show plots

printflag = 0
plotextraname = 'juanamlan_'
extralegend = ''
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
M1.updateparameter('plot_temperature',['300'])
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
M1.plotcharacterization(1,'Vth','Leff',extralegend)
if printflag==1:
  plt.savefig(plotextraname+'vthlin.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'b')
M1.plotcharacterization(2,'Vth','Leff',extralegend)
if printflag==1:
  plt.savefig(plotextraname+'vthsat.png', bbox_inches='tight')

M1.updateparameter('symbol' , 's')
M1.updateparameter('color' , 'b')
M1.plotcharacterization(3,'DIBL','Leff',extralegend)
axes = plt.gca()
axes.set_xlim([0.03,1])
if printflag==1:
  plt.savefig(plotextraname+'dibl.png', bbox_inches='tight')


M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('color' , 'b')
M1.updateparameter('symbol' , 's')
M1.plotcharacterization(4,'SS','Leff',extralegend)
axes = plt.gca()
axes.set_ylim([0,2])
if printflag==1:
  plt.savefig(plotextraname+'sslin.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('color' , 'b')
M1.updateparameter('symbol' , 'o')
M1.plotcharacterization(5,'SS','Leff',extralegend)
axes = plt.gca()
axes.set_ylim([0,2])
if printflag==1:
  plt.savefig(plotextraname+'sssat.png', bbox_inches='tight')


M1.updateparameter('plot_characterization_fit' , 1)
M1.updateparameter('plot_characterization_vdref' , 0.05)
M1.updateparameter('xlogflag' , 0)
M1.updateparameter('color' , 'b')
M1.plotcharacterization(6,'Ron','Leff',extralegend)
axes = plt.gca()
axes.set_ylim([0,1800])
#axes.set_xlim([0,0.5])
if printflag==1:
  plt.savefig(plotextraname+'ron.png', bbox_inches='tight')

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'b')
M1.plotcharacterization(7,'Gmmax','Leff',extralegend)
axes = plt.gca()
axes.set_ylim([0,0.005])
if printflag==1:
  plt.savefig(plotextraname+'gmax.png', bbox_inches='tight')
########################################show plots

M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('ylogflag' , 1)
M1.updateparameter('xlogflag' , 0)
M1.updateparameter('color' , 'b')
M1.plotcharacterization(8,'Ioff','Ion',extralegend)
if printflag==1:
  plt.savefig(plotextraname+'IonIoff_juanamlan.png', bbox_inches='tight')
  
  

M1.updateparameter('plot_characterization_save' , 1)
M1.updateparameter('plot_characterization_file_out' , '/home/juan/research/intern/SC_III-V/'+'Iofffloor_vs_Ronmax_SAT_300K.txt')
  
M1.updateparameter('plot_characterization_fit' , 0)
M1.updateparameter('plot_characterization_vdref' , 0.5)
M1.updateparameter('ylogflag' , 1)
M1.updateparameter('xlogflag' , 1)
M1.updateparameter('color' , 'b')
M1.updateparameter('legend_loc' , 'upper right')

M1.plotcharacterization(9,'Iofffloor','Ronmax',extralegend)  

plt.show() 




