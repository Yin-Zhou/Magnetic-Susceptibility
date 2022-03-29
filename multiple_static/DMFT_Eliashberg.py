from scipy.optimize import root, diagbroyden, excitingmixing
import numpy as np
from numpy.linalg import eig
import cmath


# consider we have one q point
# read T, CHI_S, CHI_C, U, mu, nr, sigi(2*Nc)
# kx, ky, kz, Nc
# ham(nr), tran(nr,3)

f_in = open('DMFT_chi_static_4.in', 'r')
lines = f_in.readlines()

numq = (lines[0]).split()[0]
T = float((lines[2]).split()[0])
Nc = int((lines[3]).split()[0])
print(Nc)
mu = float((lines[4]).split()[0])
[kx,ky,kz] = (lines[5]).split()
kx = int(kx)
ky = int(ky)
kz = int(kz)
print(kx,ky,kz)
nr = int((lines[6]).split()[0])
# print(nr)
norb = int((lines[7]).split()[0])

trans = lines[9].split()
# print(trans)
tran = np.zeros([nr,3])
for r in range(nr):
	for i in range(3):
		tran[r][i] = int(trans[3*r+i])
hams = lines[10].split()
# ['(', 'num', ',', 'num', ')']
ham = np.empty(nr,dtype=complex)
for r in range(nr):
	ham[r] = complex(float(hams[5*r+1]), float(hams[5*r+3]))
sigis = lines[11].split()
# print(sigis[:10])
sigi = np.empty(2*Nc,dtype=complex)
for i in range(2*Nc):
	sigi[i] = complex(float(sigis[5*i+1]), float(sigis[5*i+3]))
# print(lines[12].split())
U_real = float((lines[12].split())[1])
U_imag = float((lines[11].split())[3])
U = complex(U_real, U_imag)

f_in.close()

f_out = open('DMFT_chi_static_4.out', 'r')

lines = f_out.readlines()
CHI_S = float((lines[2].split())[-1])
CHI_C = float((lines[3].split())[-1])

f_out.close()

# T = 0.01
# mu = 1.0
# nr = 1
# ham = np.ones(nr)
# tran = np.ones((nr,3))
# Nc = 1
# sigi = np.ones(2*Nc)
# CHI_S = 1.
# CHI_C = 1.
# U = 1.0

# the effective interaction
Veff = U+(3./2.)*U*U*CHI_S-(1./2.)*U*U*CHI_C

# Green function (recalculate?)
# dimension of G: kx*ky*kz*2*Nc
# {tran()}
# kx=10
# ky=10
# kz=10
G = np.ones((kx*ky*kz*2*Nc,1),dtype=complex)

# ham nr, tran (nr,3)
ommesh = np.array([(2*i+1)*np.pi*T for i in range(2*Nc)])

# set up the Green function
for ikx in range(kx):
	for iky in range(ky):
		for ikz in range(kz):
			hk = 0.0
			for r in range(nr):
				hk += ham[r]*np.exp(complex(0.,1.)*(np.dot([kx,ky,kz], tran[r,:])))
			for n in range(2*Nc):
				index = ikx*ky*kz*2*Nc + iky*kz*2*Nc + ikz*2*Nc + n
				temp = complex(0.,1.)*ommesh[n] + mu - hk - sigi[n]
				G[index,0] = 1/temp

# Eliashberg equation in matrix form
# consider what Veff is, Veff(k-k_prime), but Veff(CHI_S,CHI_C,U)
A = np.empty((kx*ky*kz*2*Nc,kx*ky*kz*2*Nc),dtype=complex)
for ikx in range(kx):
	for iky in range(ky):
		for ikz in range(kz):
			for n in range(2*Nc):
				index = ikx*ky*kz*2*Nc + iky*kz*2*Nc + ikz*2*Nc + n
				A[:,index] = Veff*G[index,0]*G[index,0]
A = (-T/(kx*ky*kz))*A

identity = np.eye(A.shape[0], A.shape[1],dtype=complex)


# calculate the Eliashberg equation as a function of delta_k
def Eliashberg(delta_k):
	# delta_k should be a vector of dimension (kx*ky*kz*2*Nc,1)
	output = np.dot((A-identity), delta_k)
	return output


# the derivative of the function
# raise problem as we have to use the matrix here?
def Eliashberg_prime(delta_k):
	return (A-identity)


# calculate the largest eigenvalue
evalues, evectors = eig(A)
evalue = np.amax(evalues)
print(evalue)

# print(A.shape)
# initial guess
x0 = np.ones(kx*ky*kz*2*Nc,dtype=complex)
# find delta_k
answer = excitingmixing(Eliashberg, x0)

# print(answer.shape)
print(answer[0])









