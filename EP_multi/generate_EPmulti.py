import numpy as np
import ast

T = 0.01
###########
# mu = 8.0
###########
nu = 12
delta = 1.0e-4
norb = 5
eps_acustic = 0.00062 # in eV 
strength = 0.0

numomega = 20001
omega_lower = 0.0
omega_upper = 2.0

# numk = 6
# numq = 32
numnu = 12

# enk = [8.4735]
# enkq = [8.4735,8.4914,9.5449,10.1697,9.5449,10.1697,12.5216,12.5562]

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

numk = 125
numq = 125

enk = np.zeros((int(numk),int(norb)))
enkq = np.zeros((int(numk),int(numq),int(norb),int(norb)))

phonon_frequency = np.zeros((int(numq),int(numnu)))

# electron-phonon matrix
Gv = np.zeros((int(numk),int(numq),int(numnu),int(norb),int(norb)))

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
		start = starts[iq] + 2

		# start at the current k point
		start += 1

		# qpoints.append(lines[start].split()[-3:])
		# # start at k points
		# start += 1

		for ik in range(int(numk)):
			# if iq == 0:
			# 	kpoints.append(lines[start].split()[-3:])

			if iq == 0:
				for iorb in range(int(norb)):
					enk[ik,iorb] = lines[start+3+iorb*numnu*norb].split()[-4]

			for korb in range(int(norb)):
				for qorb in range(int(norb)):
					enkq[ik,iq,korb,qorb] = lines[start+3+qorb*numnu+korb*norb*numnu].split()[-3]

					for inu in range(int(numnu)):
						Gv[ik,iq,inu,korb,qorb] = float(lines[start+3+inu+qorb*numnu+korb*norb*numnu].split()[-1])/1000.0 #convert to eV

			if ik == 0:
				for inu in range(int(numnu)):
					phonon_frequency[iq,inu] = float(lines[start+3+inu].split()[-2])/1000.0 #convert to eV


			# for inu in range(int(numnu)):
				# if inv == 0:
				# # 	print(lines[start+3])
				# # 	print(lines[start+3].split())
				# 	enk.append(lines[start+3].split()[-4])
				# 	enkq.append(lines[start+3].split()[-3])
				# print(lines[start])
				# if ik == 0:
				# 	# print(lines[start+3+inv].split())
				# 	phonon_frequency[iq,inv] = float(lines[start+3+inv].split()[-2])/1000.0 #convert to eV
				# print(lines[start+3+inv].split()[-1])
				# Gv[ik,iq,inu] = float((lines[start+3+inu].split())[-1])/1000.0 #convert to eV

   #          # start at the next kpoint
			# start += 5 + int(numnu)
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
EP_INPUT.write('%s\n'%numnu)

for ik in range(int(numk)):
	for iorb in range(int(norb)):
		EP_INPUT.write('%s '%str(enk[ik][iorb]))
EP_INPUT.write('\n')

for ik in range(int(numk)):
	for iq in range(int(numq)):
		for korb in range(int(norb)):
			for qorb in range(int(norb)):
				EP_INPUT.write('%s '%str(enkq[ik][iq][korb][qorb]))
EP_INPUT.write('\n')

for ik in range(int(26)):
	for i in range(3):
		EP_INPUT.write('%s '%str(kpoints[ik][i]))
EP_INPUT.write('\n')

for iq in range(int(9126)):
	for i in range(3):
		EP_INPUT.write('%s '%str(qpoints[iq][i]))
EP_INPUT.write('\n')

for iq in range(int(numq)):
	for inv in range(int(numnu)):
		EP_INPUT.write('%s '%phonon_frequency[iq][inv])
EP_INPUT.write('\n')

for ik in range(int(numk)):
	for iq in range(int(numq)):
		for inv in range(int(numnu)):
			for iorb in range(int(norb)):
				for jorb in range(int(norb)):
					EP_INPUT.write('%s '%Gv[ik][iq][inv][iorb][jorb])
EP_INPUT.write('\n')

for iq in range(int(numq)):
	EP_INPUT.write('%s '%str(weights[iq]))
EP_INPUT.write('\n')

EP_INPUT.close()