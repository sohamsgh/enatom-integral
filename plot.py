import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

fig=plt.figure(figsize=(4,4), dpi= 300)    
ax = fig.add_subplot(111)    
mpl.rc('font', family='sans-serif')    
mpl.rc('font', serif='Helvetica Neue')    
mpl.rc('text', usetex='false') 
mpl.rcParams['axes.linewidth'] = 2.0
a = np.loadtxt("dos_free.dat", skiprows=1)

title = "dos_free-2.png"
ax.set_xlabel("Energy $\epsilon$ ($\\frac{4\pi^{2}}{a^{2}}$)", fontsize=12)
ax.set_ylabel("DOS (1/$\epsilon$)", fontsize=12)
ax.scatter(a[:,0], a[:,1], s=2, color="black", label="$\epsilon = \\frac{1}{2}\left ( k_{x}^{2}+k_{y}^{2}+k_{z}^{2} \\right ) $")  
#ax.scatter(a[:,0], a[:,1], s=2, color="black", label="$\sqrt {k_{x}^{2}+k_{y}^{2}+k_{z}^{2}} $")  
#ax.scatter(a[:,0], a[:,1], s=2, color="black", label="$-2\left ( \cos(2\pi k_{x}) + cos(2\pi k_{y}) + cos(2\pi k_{z} \\right )$")  
ax.plot(a[:,0], a[:,2], color="red", label="$\\frac{2\Omega\cdot 2^{3/2}}{4\pi^{2}}\sqrt{\epsilon}$") 
#ax.set_xlim(0.0, 0.10125)
#ax.set_ylim(0.0, 35.00)
plt.legend()
plt.savefig(title, bbox_inches="tight")   


