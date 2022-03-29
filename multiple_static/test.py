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

with open('spectral.out') as f:
	lines = f.readlines()
	spectrals = [line.split()[-1] for line in lines]

distance = np.linspace(0.,0.5,6)
omega = np.linspace(-1,1,50)
X,Y = np.meshgrid(distance,omega)
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
		R[i][j] = spectrals[j*50+i]

fig, ax = plt.subplots()
cs = ax.contourf(X,Y,R,100)
ax.tick_params('both',left=True,right=True,top=True,direction='in',width=3,length=7.5,)
# ax.set_xticks(paths)
# ax.set_xticklabels(symbols)
# ax.set_yticks(paths)
# ax.set_yticklabels(symbols)
ax.set_xlabel('distance')
ax.set_ylabel(r'$\omega$')
# ax.grid()
# ax.set_title(r'$\delta = {:.2f}$'.format(infsimal))
cbar = fig.colorbar(cs,ax=ax)
plt.show()



