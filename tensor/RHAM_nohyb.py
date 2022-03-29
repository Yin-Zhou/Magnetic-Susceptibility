# program for changing RHAM to RHAM_nohyb
# norb = 4
norb = 4

f = open('RHAM', 'r')
# f_nohyb = open(RHAM_nohyb, 'w')
string_list = f.readlines()

get_zeros = string_list[2].split(',')
# print('Afte split:', nohyb)
zero_to_use = get_zeros[1]
zero_to_use_end = get_zeros[-2]
zero_to_use_start = (string_list[3].split(','))[0]

num_iter = int((len(string_list) - 1)/6)
# print(num_iter) 
# num of iterations = 1331

for i in range(num_iter-1):
	index = (2 + 6 * (i + 1))
	origin = string_list[index].split(',')

	# change first row to nohyb
	for i in range(1, norb):
		if i != (norb - 1):
			origin[i] = zero_to_use
		else:
			origin[i] = zero_to_use_end

	string_list[index] = (',').join(origin)

f.close()
# f_nohyb.close()

