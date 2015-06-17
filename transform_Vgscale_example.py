#script example to transform scale of data, this targets to transform I-V data with different Vg scale
#Juan Pablo Duarte

rootfolder = '/home/juan/research'

#indicate path for folders containing required classes
import sys
import matplotlib.pyplot as plt
import measurementtests as mt

pathfolder = '/home/juan/research/intern/data/'
pathdatafile = 'KL35LC12-S0118_16FET05Idgs@Vg51d3.txt'
pathdatafileout = 'KL35LC12-S0118_16FET05Idgs@Vg51d3SCALED.txt'

###################################creating class mt
M1 = mt.mt()

####setting up parameters for scale transformation
M1.updateparameter('var_scale_ref','Vg`') #ex: Vg, then data is y=f(Vg)
M1.updateparameter('var_dependent',['Ig','Id','Is']) #  ex: Id, Ig, Is
M1.updateparameter('var_fixed',['Vd`'])# ex: Vd

#scale reference: initial, final and delta values
M1.updateparameter('var_scale_ref_init',-5.000e-01)
M1.updateparameter('var_scale_ref_final',1.0)
M1.updateparameter('var_scale_ref_delta',0.035)
M1.convertscale(pathfolder + pathdatafile,pathfolder+pathdatafileout)

##################################plot of both results
M1.updateparameter('filetype','SAS')#update format type of data to be input
#plot
M1.updateparameter('symbol','o') 
M1.updateparameter('markersize',5)
M1.updateparameter('color','r')
M1.plotfiledata(pathfolder + pathdatafile,'Vg`','Id',1)

M1.updateparameter('symbol','o') 
M1.updateparameter('markersize',5)
M1.updateparameter('color','b')
M1.plotfiledata(pathfolder + pathdatafileout,'Vg`','Id',1)

plt.show() 



