#analyze Ioff data
#Juan Duarte
#SC III-V

from numpy import loadtxt
import matplotlib.pyplot as plt
import numpy as np

pathandfile = '/home/juan/research/intern/SC_III-V/summaryIoffRon.txt'

datalist = loadtxt(pathandfile,skiprows = 1)
devicenumber = datalist[:,0]
Leffall = datalist[:,1]

Ioff400K = datalist[:,6]
Ioff300K = datalist[:,4]
Ioff250K = datalist[:,8]
Ioff200K = datalist[:,10]
Ioff133K = datalist[:,12]

Leff = list(set(Leffall))
Leff = sorted(Leff)

Ldata = {}
for L in Leff:
  Ldata[L]={'133K':[],'200K':[],'250K':[],'300K':[],'400K':[]}

count=0 
for L in Leffall:
  if Ioff133K[count]!=0:
    Ldata[L]['133K'].append(Ioff133K[count])
  if Ioff200K[count]!=0:
    Ldata[L]['200K'].append(Ioff200K[count])    
  if Ioff250K[count]!=0:
    Ldata[L]['250K'].append(Ioff250K[count])
  if Ioff300K[count]!=0:
    Ldata[L]['300K'].append(Ioff300K[count])
  if Ioff400K[count]!=0:
    Ldata[L]['400K'].append(Ioff400K[count])
    
  count+=1
  
plt.figure(1)

symbolarray = ['o','o' ,'o' ,'o' ,'s' ,'s','s','s','>','>','>','>']
colorarray = ['r','k','b','g','r','k','b','g','r','k','b','g' ]
countall=0
for L in Leff:
  average = []
  oneovert = []
  for key in Ldata[L]:
    l = []
    l = Ldata[L][key] 
    if float(len(l))!=0:
      logave = map(np.log,l)
      ave= np.exp(((sum(logave))) / float(len(l)))
      average.append(ave)
      oneovert.append( 1/float(key[:-1]))
      
      namelegend = L 
  plt.plot(oneovert,average,symbolarray[countall],lw=8,markersize=10, color=colorarray[countall],label=namelegend)
  countall+=1
plt.legend(loc='upper right')

ax = plt.gca()
ax.legend(bbox_to_anchor=(1.15, 1.05))  
ax.set_yscale('log')  
ax.set_xlabel('1/T (1/K)')
ax.set_ylabel('Ioff Floor Sat (A/um)')
plt.savefig('IoffversisoneT.png', bbox_inches='tight')
plt.show()   
