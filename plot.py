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
import os;

filename = sys.argv[1];
A=loadtxt(filename, skiprows=6);

# read the energy
f=open(filename,'r')
#skip the first two lines
f.readline();
f.readline();
energystring = f.readline()
f.close()
El=energystring.split('\t')
E=El[2];

f = figure();
ad = f.add_subplot(111);

plot(A[:,0],A[:,1])

xlabel('blocksize');
ylabel('standard deviation');
title('blocking results');
text(.85,.1,'Energy: '+E,horizontalalignment='right',verticalalignment='center',transform= ad.transAxes);



show(block=False)





# parameter optimization


Parameters = loadtxt('../Proj1/parameters.txt');
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
plt.show(block=False)


# single particle density

r=[];


datapoints = os.path.getsize('../Proj1/positions.bin')/8;
p=open('../Proj1/positions.bin',mode='rb');
values = array.array('d');
values.read(p,datapoints)


for p in range(0,datapoints):
	r.append(values[p]);                 


h= plt.figure();
had=h.add_subplot(111);
n, bins, patches = hist(values, 200, normed=1,histtype='stepfilled')
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)



xlabel('r [a.u.]')
ylabel('density')
title('single particle density for hydrogen molecule')
plt.show(block=True);

quit()



