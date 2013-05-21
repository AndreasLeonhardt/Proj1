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

filename = sys.argv[1];

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
	print(E);
	E.append(float(El[3]));
	R.append(R_Start+R_Step*r);


g = plt.figure();
gad = g.add_subplot(111);
gp=plt.plot(R,E);
xlabel('Radius [a.u.]');
ylabel('Energy');
title('Energy dependence of distance between nuclei');
plt.show(block=False);


r = argmin(E)
A=loadtxt('{0}{1}{2}'.format(filename,str(r),'.txt'), skiprows=6);
f = plt.figure();
ad = f.add_subplot(111);

fp=plt.plot(A[:,0],A[:,1])

xlabel('blocksize');
ylabel('standard deviation');
title('blocking results');
text(.85,.1,'Energy: '+str(E[r]),horizontalalignment='right',verticalalignment='center',transform= ad.transAxes);



plt.show(block=True)

quit()
