
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
import matplotlib.colors as mcol
from matplotlib.ticker import MaxNLocator,FixedLocator,IndexLocator, MultipleLocator, FormatStrFormatter
import matplotlib as mpl



def readfile(fname):
	a=[];b=[];	 emptyline1 = True;
	fin = open(fname,'r');
	for line in fin.readlines():
		if(len(line.split()) > 0):
			a.append( [ float (x) for x in line.split() ] );
			emptyline1 = True;
		elif(emptyline1):
			b.append(np.array(a));
			a=[];
			emptyline1 = False;		
	return np.array(b)



a = readfile("nog/current.out")
b = readfile("g/current.out")


fig = plt.figure(figsize=(10,8))
ax = plt.gca();

dat1 = a

#totcharge, tottime, totcharge*1.0d0/tottime, zt, ntrap, nelec
# 		0, 				1, 								2, 						3, 		4,			 5

n = 10 - a[:,0,5];
print(n)

ax.plot(n, a[:,:,2], ls='', marker='x',ms=7, color='red');
ax.plot(n, b[:,:,2], ls='', marker='x',ms=7, color='blue');
ax.plot(n, np.sum(a[:,:,2],axis=1)/30, ls='', marker='o',ms=10, color='red',label='No Coupling');
ax.plot(n, np.sum(b[:,:,2],axis=1)/30, ls='', marker='s',ms=10, color='blue',label='g = 0.3 eV');
#ax.plot(n, a[:,:,2], ls='', marker='^',ms=10, color='red');
#ax.plot(n, b[:,:,2], ls='', marker='<',ms=10, color='yellow');


ax.set_xlabel("Net charge ($e$)",fontsize=20)
ax.set_ylabel("Average Current ($e*t_h$)",fontsize=20)

ax.set_title("Organic microcavity with 10 molecules",fontsize=20)
ax.set_xlim(-10,10)
ax.set_ylim(0,3.2)

ax.text(-4,3.05,"$\omega=\omega_0=2eV$, $E.r_{nns}=1eV$, PBC",fontsize=20)

ax.legend(loc=2,frameon=False)

plt.tight_layout();

plt.savefig('current-average.pdf', format='pdf');

plt.show()



