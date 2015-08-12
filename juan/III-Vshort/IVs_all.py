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

'''
macro = '35fet_GF2_FetPS_10'# # 
device = 'FET01' #'FET7_1x0.08_NPN_M1' # 
L = '0.1'
W = '12'
Temp = '300'
devicenumber='841'
'''

macro = '35fet_GF2_FetPS_10'# # 
device = 'FET08' #'FET7_1x0.08_NPN_M1' # 
L = '1'
W = '12'
Temp = '300'
devicenumber='841'





###################################definition of class mt
#filesdata = [Temp+'K.txt']#['100K.txt','133K.txt','200K.txt','250K.txt','300K.txt','400K.txt']
#temperatures = [Temp+'K']#['100K','133K','200K','250K','300K','400K']
filesdata = ['133K_1.txt','200K_1.txt','250K_1.txt','300K_1.txt','400K_1.txt']
temperatures = ['133K','200K','250K','300K','400K']

colorarray = [['r','b','c','k','m','g'],['b','c','k','m','g','r'],['c','k','m','g','r','b'],['k','m','g','r','b','c'],['m','g','r','b','c','k'],['g','r','b','c','k','m']]

currents = ['Id']#,'Is','Ig'
vdsarray = [0.05]
count3=-1
for vds in vdsarray:
  count3+=1
  count2=-1  
  for current in currents:
    count2+=1
    count=-1
    for filename in filesdata:
      count+=1
      M1 = mt.mt()

      ###################################input data for different devices
      M1.updateparameter('device_tags',['technology','wafer','site','macro','device','test','temperature'])

      M1.add_die_info('/home/juan/research/intern/SC_III-V/SA35FET19.txt')
      #add entire
      M1.adddata('/home/juan/research/intern/SC_III-V/ALLDATA/'+filename,temperatures[count])

      #######################################plot set up
      M1.updateparameter('plot_technology','all')
      M1.updateparameter('plot_wafer','all')
      M1.updateparameter('plot_site','all')
      M1.updateparameter('plot_macro',[macro])#width,35fet_GF2_FetPS_10
      M1.updateparameter('plot_device',[device] )#->gate length ['1.0x0.08_P']  ['FET3_5x0.1_NPN_M1']1.0x0.08_P
      #['FET7_1x0.08_NPN_M1','FET5_5x0.08_NPN_M1_CAoverPC']
      M1.updateparameter('plot_test',['Idgsx@Vg26d2','Idgsx@Vg126d2'])
      M1.updateparameter('plot_temperature','all')
      M1.updateparameter('plot_y_variables',['Id'])
      M1.updateparameter('plot_x_variable','Vg`')
      M1.updateparameter('plot_legend_names',['bias','temperature'])#'Leff','Wdes'
      M1.updateparameter('plot_title_names',['technology','wafer','site','macro'])
      M1.updateparameter('plot_parameter_limits',[ ['Vd`',[vds]] ])#, ['Vd`',[1.0]]
      M1.updateparameter('plot_W_normalization_flag',1)#
      ######################################plot run
      M1.updateparameter('legend_loc' , 'upper left')
      M1.updateparameter('plot_y_variables',[current])
      M1.updateparameter('symbol' , '-')
      M1.updateparameter('lw' , 4)
      M1.updateparameter('color' , colorarray[count][count3])
      M1.updateparameter('plot_all_together' , 1)
      M1.updateparameter('ylogflag' , 0)
      M1.plotdevices((count2)*2+1)



ax = plt.gca()
plt.title(M1.plot_macro[0]+', FET L='+L+', W='+W)
ax.set_ylabel(currents[0]+'/um (A/um)')
#plt.savefig(Temp+'K_lin_'+ M1.plot_macro[0]+'_'+device+'_'+currents[0], bbox_inches='tight')

count3=-1
for vds in vdsarray:
  count3+=1
  count2=-1  
  for current in currents:
    count2+=1
    count=-1
    for filename in filesdata:
      count+=1
      M1 = mt.mt()

      ###################################input data for different devices
      M1.updateparameter('device_tags',['technology','wafer','site','macro','device','test','temperature'])

      M1.add_die_info('/home/juan/research/intern/SC_III-V/SA35FET19.txt')
      #add entire
      M1.adddata('/home/juan/research/intern/SC_III-V/ALLDATA/'+filename,temperatures[count])

      #######################################plot set up
      M1.updateparameter('plot_technology','all')
      M1.updateparameter('plot_wafer','all')
      M1.updateparameter('plot_site','all')
      M1.updateparameter('plot_macro',[macro])#width,35fet_GF2_FetPS_10
      M1.updateparameter('plot_device',[device] )#->gate length ['1.0x0.08_P']  ['FET3_5x0.1_NPN_M1']1.0x0.08_P
      #['FET7_1x0.08_NPN_M1','FET5_5x0.08_NPN_M1_CAoverPC']
      M1.updateparameter('plot_test',['Idgsx@Vg26d2','Idgsx@Vg126d2'])
      M1.updateparameter('plot_temperature','all')
      M1.updateparameter('plot_y_variables',['Id'])
      M1.updateparameter('plot_x_variable','Vg`')
      M1.updateparameter('plot_legend_names',['bias','temperature'])#'Leff','Wdes'
      M1.updateparameter('plot_title_names',['technology','wafer','site','macro'])
      M1.updateparameter('plot_parameter_limits',[  ['Vd`',[vds]]])#, ['Vd`',[1.0]] 
      M1.updateparameter('plot_W_normalization_flag',1)#
      ######################################plot run
      M1.updateparameter('lw' , 4)
      M1.updateparameter('color' , colorarray[count][count3])
      M1.updateparameter('plot_all_together' , 1)
      M1.updateparameter('legend_loc' , 'lower right')
      M1.updateparameter('plot_y_variables',[current])
      M1.updateparameter('symbol' , '-')
      M1.updateparameter('ylogflag' , 1)
      M1.plotdevices((count2)*2+2)

ax = plt.gca()
plt.title(M1.plot_macro[0]+', FET L='+L+', W='+W)
ax.set_ylabel(currents[0]+'/um (A/um)')
#plt.savefig(Temp+'K_log_'+ M1.plot_macro[0]+'_'+device+'_'+currents[0], bbox_inches='tight')



plt.show() 




