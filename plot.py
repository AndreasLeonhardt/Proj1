# Python script for plotting the result of the calculation
# called at the end of the programm Proj1.cpp
# displays the variance over the block size, as calculated in mcInt.blocking(..) 
# and writes the energy in the plot. 
# takes the data from the result file, given in the parameter file 'parameters.cfg' by the keyword 'outputfile
# at current stage it skips the first few lines according to the number specified below.
# Any modifications in creating the result file have to take care of this.
# also the position of the energy is important (3rd line, 3rd position)
# created 19.04.2013

from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import array
import os

filename = sys.argv[1];


# Energy-R dependence plot ------------------------------------------------------------------------
rmax = int(sys.argv[2]);
R_Start = float(sys.argv[3]);
R_Step = float(sys.argv[4]);


E=[];
R=[];
for r in range(0,rmax):
	# read the energy
	f=open('{0}{1}{2}'.format(filename,r,'.txt'));
	#skip the first two lines
	f.readline();
	f.readline();
	energystring = f.readline()
	f.close()
	El=energystring.split('\t')
	E.append(float(El[3]));
	R.append(R_Start+R_Step*r);


g = plt.figure();
gad = g.add_subplot(111);
gp=plt.plot(R,E);
xlabel('Radius [a.u.]');
ylabel('Energy [a.u.]');
title('Energy dependence of distance between nuclei');
plt.show(block=False);


# Single particle density -------------------------------------------------------------------------
#open positions ( not dynamically)
r=[];
x=[];
z=[];

R0 = argmin(E)
datapoints = os.path.getsize('../Proj1/samples_R{0}_positions_0.dat'.format(str(R0)))/8;
p=open('../Proj1/samples_R{0}_positions_0.dat'.format(str(R0)),mode='rb');
values = array.array('d');
values.read(p,datapoints)


for p in range(0,datapoints/3):
	r.append(values[3*p]);                 
 	x.append(values[3*p+1]);               
 	z.append(values[3*p+2]);

h= plt.figure();
had=h.add_subplot(111);
twoDhist,xedges,yedges = histogram2d(x,z,bins=200);
twoDhist/=datapoints/3;
twoDhist.shape, xedges.shape, yedges.shape
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
plt.imshow(twoDhist, extent=extent, interpolation='nearest')
plt.colorbar()
plt.show(block=False)
xlabel('z [a.u.]')
ylabel('x [a.u.]')
title('single particle density for hydrogen molecule')
plt.show(block=False);



# Blocking plot ------------------------------------------------------------------------------------
A=loadtxt('{0}{1}{2}'.format(filename,str(R0),'.txt'), skiprows=6);
f = plt.figure();
ad = f.add_subplot(111);

fp=plt.plot(A[:,0],A[:,1])

xlabel('blocksize');
ylabel('energy standard deviation [a.u.]');
title('blocking results');
text(.85,.1,'Energy: '+str(E[R0])+'a.u.',horizontalalignment='right',verticalalignment='center',transform= ad.transAxes);

plt.show(block=False)


# parameter plot ------------------------------------------------------------------------------------

Parameters = loadtxt('../Proj1/optimization_R{0}.txt'.format(str(R0)));
k =plt.figure();
ad = k.add_subplot(111)
kp=plt.plot(Parameters[:,0])
xlabel('iteration')
ylabel('alpha')
title('parameter optimization for alpha')
plt.show(block=False)

l =plt.figure();
ad = k.add_subplot(111)
lp=plt.plot(Parameters[:,1])
xlabel('iteration')
ylabel('beta')
title('parameter optimization for beta')
plt.show(block=True)


quit()
