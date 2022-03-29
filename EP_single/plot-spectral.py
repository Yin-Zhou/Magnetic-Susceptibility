import numpy as np 
import matplotlib.pyplot as plt 

infsimal = 1e-1

# Pauli matrices
tau0 = np.eye(2,dtype=np.complex128)
tau3 = np.array([[1,0],[0,-1]],dtype=np.complex128)
tau1 = np.array([[0,1],[1,0]],dtype=np.complex128)

delta0 = 0.035 # eV
deltak = lambda kx,ky: delta0*(np.cos(kx)-np.cos(ky))/2

t1 = 0.38 # eV
t2 = 0.32*t1
t3 = 0.5*t2 

epsilonk = lambda kx,ky: -2*t1*(np.cos(kx)+np.cos(ky))+4*t2*np.cos(kx)*np.cos(ky)-2*t3*(np.cos(2*kx)+np.cos(2*ky))

green = lambda kx,ky,ome: (ome+infsimal*1j)*tau0+epsilonk(kx,ky)*tau3+deltak(kx,ky)*tau1
greensc = lambda kx,ky,ome: green(kx,ky,ome)/((ome+infsimal*1j)**2-np.sqrt(epsilonk(kx,ky)**2+deltak(kx,ky)**2))

N = 101
M = 101
omega = np.linspace(-5,5,N)
# omega = 0
kpoints = np.linspace(0,np.pi,M)
ky = np.linspace(0,np.pi,M)

# test1
# spectral = lambda kx,ome: -np.imag(1/(ome-epsilonk(kx,0)+infsimal*1j))/np.pi

# test2
spectral = lambda kx,ky,ome: -np.imag(greensc(kx,ky,ome)[0][0])/np.pi

symbols = [r'$\Gamma$', 'X']
paths = [kpoints[0],kpoints[-1]]

with open('Ksym.txt') as f:
	lines = f.readlines()
	numk = lines[0].split()[0]

with open('spectral.out') as f:
	lines = f.readlines()
	distances = [line.split()[0] for line in lines]
	omegas = [line.split()[1] for line in lines]
	spectrals = np.array([line.split()[-1] for line in lines]).astype(np.float)

num_omega = 20001
# omega = np.linspace(0,2,num_omega)
omega = []
for i in range(num_omega):
	omega.append(omegas[i])
omega = np.array(omega).astype(np.float)
distance = []
for i in range(int(numk)):
	distance.append(distances[i*num_omega])
distance = np.array(distance).astype(np.float)
paths = [distance[0],distance[-1]]
# print(paths)
X,Y = np.meshgrid(distance,omega)
# print(X[:10])
# print(Y[:10])
# X,Y = np.meshgrid(kpoints,ky)
# print(len(spectrals))
# print(len(X))
# print(X)
# print(len(X[0]))
# print(X.shape)
R = np.empty(X.shape)
for i in range(len(X)):
	for j in range(len(X[0])):
		# xcoord = X[i][j]
		# ycoord = Y[i][j]
		# R[i][j] = spectral(xcoord,ycoord,omega)
		R[i][j] = spectrals[j*num_omega+i]

fig, ax = plt.subplots()
cs = ax.contourf(X,Y,R,100)
ax.tick_params('both',left=True,right=True,top=True,direction='in',width=3,length=7.5,)
ax.set_xticks(paths)
ax.set_xticklabels([r'$\Gamma$',r'$X$'])
# ax.set_yticks(paths)
# ax.set_yticklabels(symbols)
ax.set_xlabel('distance')
ax.set_ylabel(r'$\omega$')
ax.set_ylim(-2,-0.5)
# ax.grid()
# ax.set_title(r'$\delta = {:.2f}$'.format(infsimal))
cbar = fig.colorbar(cs,ax=ax)
# plt.locator_params(axis='x', nbins=5)
# fig.savefig('spectral.pdf',bbox_inches='tight',dpi=300)
# plt.close('all')
plt.show()



