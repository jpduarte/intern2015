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

#FET 1P, 4
#large dependance of Vd on Ig
macro = 'COLADA_JTK_Device1P'# # 
device = '0.16x0.20_P' #'FET7_1x0.08_NPN_M1' # 
L = '0.2'
W = '0.16'
Temp = '300'
devicenumber='4'

'''
#FET 2P, 5
macro = 'COLADA_JTK_Device2P'# # 
device = '1.0x0.08_P' #'FET7_1x0.08_NPN_M1' # 
L = '0.08'
W = '1'
Temp = '300'
devicenumber='5'
'''

'''
#FET 2P, 7
macro = 'COLADA_JTK_Device2P'# # 
device = '1.0x0.14_P' 
L = '0.14'
W = '1'
Temp = '300'
devicenumber='7'
'''

'''
#BJT NM1,7
#little Ig dependance on Vd
macro = 'NPN_M1'# # 
device ='FET7_1x0.08_NPN_M1' # 
L = '0.08'
W = '1'
Temp = '300'
devicenumber='5'
'''
'''
#BJT NM1,5
macro = 'NPN_M1'# # 
device ='FET5_5x0.08_NPN_M1_CAoverPC' # 
L = '0.08'
W = '5'
Temp = '300'
devicenumber='5'
'''
'''
macro = 'PNP_M1'
device ='FET7_1x0.08_PNP_M1'
L = '0.08'
W = '1'
Temp = '300'
devicenumber='7'
'''



###################################definition of class mt
#filesdata = [Temp+'K.txt']#['100K.txt','133K.txt','200K.txt','250K.txt','300K.txt','400K.txt']
#temperatures = [Temp+'K']#['100K','133K','200K','250K','300K','400K']
filesdata = ['100K.txt','133K.txt','200K.txt','250K.txt','300K.txt','400K.txt']
temperatures = ['100K','133K','200K','250K','300K','400K']

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

      M1.add_die_info('/home/juan/research/intern/III-V-ART/ARTinfotext.txt')
      #add entire
      M1.adddata('/home/juan/research/intern/III-V-ART/III-V-ART_final/'+filename,temperatures[count])

      #######################################plot set up
      M1.updateparameter('plot_technology','all')
      M1.updateparameter('plot_wafer','all')
      M1.updateparameter('plot_site','all')
      M1.updateparameter('plot_macro',[macro])#width,35fet_GF2_FetPS_10
      M1.updateparameter('plot_device',[device] )#->gate length ['1.0x0.08_P']  ['FET3_5x0.1_NPN_M1']1.0x0.08_P
      #['FET7_1x0.08_NPN_M1','FET5_5x0.08_NPN_M1_CAoverPC']
      M1.updateparameter('plot_test',['IdVgNFET_floating'])
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
      M1.plotdevices((count2)*2+1)

      M1.updateparameter('legend_loc' , 'lower right')
      M1.updateparameter('plot_y_variables',[current])
      M1.updateparameter('symbol' , '-')
      M1.updateparameter('ylogflag' , 1)
      #M1.plotdevices((count2)*2+2)



ax = plt.gca()
plt.title(M1.plot_macro[0]+', FET L='+L+', W='+W)
ax.set_ylabel(currents[0]+'/um (A/um)')
#plt.savefig(Temp+'K_lin_'+ M1.plot_macro[0]+'_'+devicenumber+'_'+currents[0], bbox_inches='tight')
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

      M1.add_die_info('/home/juan/research/intern/III-V-ART/ARTinfotext.txt')
      #add entire
      M1.adddata('/home/juan/research/intern/III-V-ART/III-V-ART_final/'+filename,temperatures[count])

      #######################################plot set up
      M1.updateparameter('plot_technology','all')
      M1.updateparameter('plot_wafer','all')
      M1.updateparameter('plot_site','all')
      M1.updateparameter('plot_macro',[macro])#width,35fet_GF2_FetPS_10
      M1.updateparameter('plot_device',[device] )#->gate length ['1.0x0.08_P']  ['FET3_5x0.1_NPN_M1']1.0x0.08_P
      #['FET7_1x0.08_NPN_M1','FET5_5x0.08_NPN_M1_CAoverPC']
      M1.updateparameter('plot_test',['IdVgNFET_floating'])
      M1.updateparameter('plot_temperature','all')
      M1.updateparameter('plot_y_variables',['Id'])
      M1.updateparameter('plot_x_variable','Vg`')
      M1.updateparameter('plot_legend_names',['bias','temperature'])#'Leff','Wdes'
      M1.updateparameter('plot_title_names',['technology','wafer','site','macro'])
      M1.updateparameter('plot_parameter_limits',[  ['Vd`',[vds]]])#, ['Vd`',[1.0]] 
      M1.updateparameter('plot_W_normalization_flag',1)#
      ######################################plot run
      M1.updateparameter('legend_loc' , 'upper left')
      M1.updateparameter('plot_y_variables',[current])
      M1.updateparameter('symbol' , '-')
      M1.updateparameter('lw' , 4)
      M1.updateparameter('color' , colorarray[count][count3])
      M1.updateparameter('plot_all_together' , 1)
      #M1.plotdevices((count2)*2+1)

      M1.updateparameter('legend_loc' , 'lower right')
      M1.updateparameter('plot_y_variables',[current])
      M1.updateparameter('symbol' , '-')
      M1.updateparameter('ylogflag' , 1)
      M1.plotdevices((count2)*2+2)

ax = plt.gca()
plt.title(M1.plot_macro[0]+', FET L='+L+', W='+W)
ax.set_ylabel(currents[0]+'/um (A/um)')
#plt.savefig(Temp+'K_log'+ M1.plot_macro[0]+'_'+devicenumber+'_'+currents[0], bbox_inches='tight')



plt.show() 




