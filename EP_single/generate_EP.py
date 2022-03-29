import numpy as np
import ast

T = 0.01
###########
# mu = 8.0
###########
nu = 12
delta = 1.0e-4
norb = 1
eps_acustic = 0.00062 # in eV 
strength = 0.0

numomega = 20001
omega_lower = -2.0
omega_upper = 0.0

# numk = 6
# numq = 32
numnv = 12

# enk = [8.4735]
# enkq = [8.4735,8.4914,9.5449,10.1697,9.5449,10.1697,12.5216,12.5562]
enk = []
enkq = []

# kpoints = [[0.,0.,0.]]
# qpoints = [[0.,0.,0.],[0.,0.,-0.5],[0.,-0.5,0.],[0.,-0.5,-0.5],[-0.5,0.,0.],[-0.5,0.,-0.5],\
# [-0.5,-0.5,0.],[-0.5,-0.5,-0.5]]
kpoints = []
qpoints = []

# weights = [0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]
weights = []

with open('Ksym.txt') as f:
	lines = f.readlines()
	numk = lines[0].split()[0]
	kpoints = [line.split()[-3:] for line in lines[1:]]

with open('Q.txt') as f:
	lines = f.readlines()
	numq = (lines[0].split())[0]
	qpoints = [line.split()[:3] for line in lines[1:]]
	weights = [(line.split()[-1])[:-2] for line in lines[1:]]

phonon_frequency = np.zeros((int(numq),int(numnv)))

# electron-phonon matrix
Gv = np.zeros((int(numk),int(numq),int(numnv)))

lookup = 'Electron-phonon vertex |g| (meV)'

find_mu = 'Fermi energy is calculated from the fine k-mesh'

with open('epw.out') as f:
	lines = []
	starts = []
	for num,line in enumerate(f,1):
	    lines.append(line)
	    if lookup in line:
	    	# should be a list of length numq
	    	starts.append(num)
	    if find_mu in line:
	    	position_mu = num
	# input mu
	mu = lines[position_mu-1].split()[-2]
	# print(mu)

	# # start at q points
	# start = start + 2

	# print(len(starts))

	for iq in range(int(numq)):
        
        # start at the current q point
		start = starts[iq] + 1

		# start at the current k point
		start += 1

		# qpoints.append(lines[start].split()[-3:])
		# # start at k points
		# start += 1

		for ik in range(int(numk)):
			# if iq == 0:
			# 	kpoints.append(lines[start].split()[-3:])

			for inv in range(int(numnv)):
				if inv == 0:
				# 	print(lines[start+3])
				# 	print(lines[start+3].split())
					enk.append(lines[start+3].split()[-4])
					enkq.append(lines[start+3].split()[-3])
				# print(lines[start])
				if ik == 0:
					# print(lines[start+3+inv].split())
					phonon_frequency[iq,inv] = float(lines[start+3+inv].split()[-2])/1000.0 #convert to eV
				# print(lines[start+3+inv].split()[-1])
				Gv[ik,iq,inv] = float((lines[start+3+inv].split())[-1])/1000.0 #convert to eV

            # start at the next kpoint
			start += 5 + int(numnv)
		# #start at the next qpoint
		# start += 2

EP_INPUT = open('EP.in', 'w')

EP_INPUT.write('%s\n'%T)
EP_INPUT.write('%s\n'%mu)
EP_INPUT.write('%s\n'%delta)
EP_INPUT.write('%s\n'%nu)
EP_INPUT.write('%s\n'%norb)
EP_INPUT.write('%s\n'%eps_acustic)
EP_INPUT.write('%s\n'%strength)

EP_INPUT.write('%s\n'%numk)
EP_INPUT.write('%s %s %s'%(numomega,omega_lower,omega_upper))
EP_INPUT.write('\n')

EP_INPUT.write('%s\n'%numq)
EP_INPUT.write('%s\n'%numnv)

for ik in range(int(numk)):
	EP_INPUT.write('%s '%str(enk[ik]))
EP_INPUT.write('\n')

for ik in range(int(numk)):
	for iq in range(int(numq)):
		EP_INPUT.write('%s '%str(enkq[int(numk)*iq+ik]))
EP_INPUT.write('\n')

for ik in range(int(numk)):
	for i in range(3):
		EP_INPUT.write('%s '%str(kpoints[ik][i]))
EP_INPUT.write('\n')

for iq in range(int(numq)):
	for i in range(3):
		EP_INPUT.write('%s '%str(qpoints[iq][i]))
EP_INPUT.write('\n')

for iq in range(int(numq)):
	for inv in range(int(numnv)):
		EP_INPUT.write('%s '%phonon_frequency[iq][inv])
EP_INPUT.write('\n')

for ik in range(int(numk)):
	for iq in range(int(numq)):
		for inv in range(int(numnv)):
			EP_INPUT.write('%s '%Gv[ik][iq][inv])
EP_INPUT.write('\n')

for iq in range(int(numq)):
	EP_INPUT.write('%s '%str(weights[iq]))
EP_INPUT.write('\n')

EP_INPUT.close()
