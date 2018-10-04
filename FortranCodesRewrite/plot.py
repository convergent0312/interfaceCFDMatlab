axstring = 'auto'

import matplotlib.pyplot as plt
import numpy as np

#importing data
datarho = np.genfromtxt('rdata.dat',delimiter=None)
datav = np.genfromtxt('vdata.dat',delimiter=None)
datap = np.genfromtxt('prdata.dat',delimiter=None)
dataphi = np.genfromtxt('phidata.dat',delimiter=None)
datagamma = np.genfromtxt('gammadata.dat',delimiter=None)

#declaring data
#X
x = datarho[:,0]

#RHO
rho_int = datarho[:,1]
rho = datarho[:,2]

#V
vint = datav[:,1]
v = datav[:,2]


#P
pint = datap[:,1]
p = datap[:,2]

#PHI

phi_int = dataphi[:,1]
phi = dataphi[:,2]

#GAMMA

gamma_int = datagamma[:,1]
gamma = datagamma[:,2]


#setting windows
#manager = plt.get_current_fig_manager()
#manager.resize(*manager.window.maxsize())

plt.clf()


plt.subplot(231)
plt.plot(x,rho_int,'-xr')
plt.hold
plt.plot(x,rho,'-ob')
plt.axis(axstring)
plt.title('RHO')

plt.subplot(232)
plt.plot(x,phi_int,'-xr')
plt.hold
plt.plot(x,phi,'-ob')
plt.title('PHI')



plt.subplot(233)
plt.plot(x,gamma_int,'-xr')
plt.hold
plt.plot(x,gamma,'-ob')
plt.title('GAMMA')

plt.subplot(234)
plt.plot(x,pint,'-xr')
plt.hold
plt.plot(x,p,'-ob')
plt.title('P')


plt.subplot(236)
plt.plot(x,vint, '-xr')
plt.hold
plt.plot(x,v,'-ob')
plt.title('V')




plt.savefig('plot.jpg')

