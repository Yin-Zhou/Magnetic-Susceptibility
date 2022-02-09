import copy

num_insertion = 5
f_origin = open('8DMFT_chi_static.in', 'r')

if f_origin.mode == 'r':
	lines = f_origin.readlines()

	numq = lines[0]
	q = lines[1]
	T = lines[2]
	Nc = lines[3]
	mu = lines[4]
	nk = lines[5]
	nr = lines[6]
	norb = lines[7]
	ncor_orb = lines[8]
	tran = lines[9]
	ham = lines[10]
	sigi = lines[11]

	q = q.split()
	qs = []
	for i in range(0,len(q),3):
		qs.append([float(q[i]), float(q[i+1]), float(q[i+2])])

	new_qs = []
	for i in range(0,len(qs)-1):
		q_prev = qs[i]
		q_next = qs[i+1]
		gap1 = (q_next[0] - q_prev[0])/(num_insertion+1)
		gap2 = (q_next[1] - q_prev[1])/(num_insertion+1)
		gap3 = (q_next[2] - q_prev[2])/(num_insertion+1)
		new_qs.append(q_prev)
		for j in range(num_insertion):
			q_toappend = [q_prev[0]+gap1,q_prev[1]+gap2,q_prev[2]+gap3]
			new_qs.append(q_toappend)
			q_prev = q_toappend
	new_qs.append(qs[-1])
	numq = int(len(new_qs))

	f = open('DMFT_chi_static.in','w+')
	f.write(str(numq)+'\n')
	for i in range(len(new_qs)):
		for j in range(len(new_qs[0])):
			if (i == len(new_qs)-1) and (j == len(new_qs[0])-1):
				f.write(str(new_qs[i][j]))
			else:
				f.write(str(new_qs[i][j])+' ')
	f.write('\n'+T)
	f.write(Nc)
	f.write(mu)
	f.write(nk)
	f.write(nr)
	f.write(norb)
	f.write(ncor_orb)
	f.write(tran)
	f.write(ham)
	f.write(sigi)
	f_origin.close()
	f.close()

else:
	print('file not in the read mode!')








