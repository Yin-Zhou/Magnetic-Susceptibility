import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


with open('self-energy.out') as f:
	lines = f.readlines()
	distances = np.array([line.split()[0] for line in lines]).astype(np.float)
	E_minus_mu = np.array([line.split()[1] for line in lines]).astype(np.float)
	imag_sigi = 1000*np.array([line.split()[-1] for line in lines]).astype(np.float) # convert to meV

points = np.array([distances,E_minus_mu]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)

fig,ax = plt.subplots(1,1,sharex=True,sharey=True)

norm = plt.Normalize(imag_sigi.min(),imag_sigi.max())
lc = LineCollection(segments,cmap='viridis',norm=norm)
lc.set_array(imag_sigi)
lc.set_linewidth(2)
line = ax.add_collection(lc)
fig.colorbar(line,ax=ax)

paths = [distances[0],distances[-1]]
symbols = [r'$\Gamma$', r'$X$']
ax.set_title('some material some temperature', fontsize=16)
ax.set_xlabel('distance',fontsize=16)
ax.set_ylabel(r'$E-E_{F} (eV)$',fontsize=16)
ax.tick_params('both',left=True,right=True,top=True,direction='in',width=3,length=7.5,)
ax.set_xticks(paths)
ax.set_xticklabels(symbols)
ax.set_xlim(distances.min(),distances.max())
# ax.set_ylim(E_minus_mu.min()-1e-1,E_minus_mu.max()+1e-1)
ax.set_ylim(-2,-0.5)
plt.show()








