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



show()

