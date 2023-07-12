# This is an experimental PoC for gpu acceleration of algorithm 2
# This was run using a single AMD GPU using WSL2
# The last time this was executed it was stable but with a different development environment there are no guarentees

import numpy as np
import pyopencl as cl
import itertools
import time
import layout

def run_ocl_kernal(queue, kernal, global_size, input_tuples, output_tuples, local_size = (32,)):
	for array, buffer in input_tuples:
		cl.enqueue_copy(queue, src=array, dest=buffer)

	kernal_arguments = [buffer for (_, buffer) in input_tuples]
	kernal_arguments += [buffer for (_, buffer) in output_tuples]

	kernal(queue, global_size, local_size, *kernal_arguments)

	for (arr, buffer) in output_tuples:
		cl.enqueue_copy(queue, src=buffer, dest=arr)

	queue.finish()


def check_sum_results(a, b, c):
	c_ref = a + b
	err = np.abs(c - c_ref)
	if ((err.sum() > 0.0).any()):
		print("err")
	else:
		print("succ")

def create_input_memory(context, input_arrays):
	return [(array, cl.Buffer(ctx, flags=cl.mem_flags.READ_ONLY, size=array.nbytes)) for array in input_arrays]

def create_output_memory(context, output_arrays):
	return [(array, cl.Buffer(ctx, flags=cl.mem_flags.WRITE_ONLY, size=array.nbytes)) for array in output_arrays]

def displayInfo():
	device = cl.get_platforms()[0]
	

ctx = cl.create_some_context(0)

#0.155
#0.171
#0.503
#4.328
#
#

distance = 6
rotated = True
#distance 8 max for 64bit, higher distances will require different data structure

if rotated:
	qubits = distance**2
else:
	qubits = distance**2 + (distance-1)**2
stabiliser_rawX, stabiliser_rawZ = layout.get_stabs(distance, rotated, True)
stabiliser_rawXnp = []
for stab in stabiliser_rawX:
	stabiliser_rawXnp.append(np.zeros(qubits).view(np.uint64))
	for bit in stab:
		if bit != -1:
			stabiliser_rawXnp[-1][bit] = 1

stabilisersX = len(stabiliser_rawX)
stabilisersZ = len(stabiliser_rawZ)
stabiliser = np.array(stabiliser_rawZ).astype(int)
num_stabilisers = len(stabiliser)
stabdistil = []
stabdirectory = [set() for i in range(stabilisersZ)]
for i in range(qubits):
	stabdistil.append(np.count_nonzero(stabiliser==i, axis=1))

for i in range(len(stabdistil)):
	for qubit in range(len(stabdistil[i])):
		if stabdistil[i][qubit] == 1:
			stabdirectory[qubit].add(i)



trajectory = int("0101",2)

t = np.unpackbits(np.array([int(''.join(list(reversed(bin(trajectory)[2:].zfill(stabilisersZ)[-stabilisersZ:]))),2)], dtype='>i8').view(np.uint8))[-stabilisersZ:]

tx_bits = np.unpackbits(np.array([int(''.join(list(reversed(bin(trajectory)[2:].zfill(stabilisersZ+stabilisersX)))),2)], dtype='>i8').view(np.uint8))[-stabilisersX:]
trajectories = [0]

cont_bits = []
combinations = []

comb_lengths = []
used = set()



for i in range(stabilisersZ):
	combinations.append([[],[]])
	firstbit = t[i]
	stabs = set(stabdirectory[i])
	stabs.difference_update(used)
	
	for c in range(len(stabs)+1):
		combinations[i][c%2] += list(itertools.combinations(stabs, c))



	cb = [stab for stab in stabdirectory[i] if stab in used]
	cont_bits.append(cb)
	comb_lengths.append(len(stabs))

	bits_contributing_to_c = (firstbit, cb)

	#####################trajectories, bits_contributing, combinations

	# print("bit:", bits_contributing_to_c[0])
	# print("bit index to check:", bits_contributing_to_c[1])
	# print("combs0", combinations[i][0])
	# print("combs1", combinations[i][1])

	for q in stabs:
		used.add(q)

bpcache = {}
def build_program(a, b):
	key = str(a)+"|"+str(b)
	if key in bpcache:
		program_source, expCombs = bpcache[key]
		return program_source, expCombs
	program_source = """

		constant int comb_mask_mask[5] = {0, 1, 2, 4, 8};
		constant int comb_mask[5][2][8] = {
			{{-1, -1, -1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1, -1, -1}},
			{{0, -1, -1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, -1, -1, -1}},
			{{0, 2, -1, -1, -1, -1, -1, -1},{1, 1, -1, -1, -1, -1, -1, -1}},
			{{0, 2, 2, 2, -1, -1, -1, -1},{1, 1, 1, 3, -1, -1, -1, -1}},
			{{0, 2, 2, 2, 2, 2, 2, 4},{1, 1, 1, 1, 3, 3, 3, 3}}
		};
	"""




	max_combs = 0
	program_source+=("	constant int comb["+str(len(combinations))+"][2][8][4] = {\n") 
	first0 = True
	for i in range(len(combinations)):
		# for p in combinations[i]:
		# 	print("c",p)
		pad = np.zeros([2, 8, 4]).view(np.uint64)
		for j in range(len(combinations[i])):
			for k in range(len(combinations[i][j])):
				for l in range(len(combinations[i][j][k])):
					pad[j][k][l] = combinations[i][j][k][l]
		# for p in pad:
		# 	print("p",p)
		first1 = True
		if first0:
			first0 = False
			program_source += "\n{"
		else:
			program_source += ",\n{"
		for i in pad:
			first2 = True
			if first1:
				first1 = False
				program_source += "\n{"
			else:
				program_source += ",\n{"
			for j in i:
				first3 = True
				if first2:
					first2 = False
					program_source += "\n	{"
				else:
					program_source += ",\n	{"
				program_source += ", ".join([str(k) for k in list(j)])
				program_source += "}"
			program_source += "}"


		program_source += "}"
	program_source += "};\n"




	cont_bits_lengths = []
	program_source += ("	constant int cont_bits["+str(len(cont_bits))+"][2] = {")
	first = True
	for bit in range(len(cont_bits)):
		# print(cont_bits[bit])
		if first:
			first = False
			program_source += "\n"
		else:
			program_source += ",\n"
		pad = [-1, -1]
		for i in range(len(cont_bits[bit])):
			pad[i] = cont_bits[bit][i]


		program_source += "{" + ', '.join([str(l) for l in pad]) + "}"
		cont_bits_lengths.append(len(cont_bits[bit]))

	program_source += ("};\n")

	program_source+=("	constant int cont_bits_lengths["+str(len(cont_bits_lengths))+"] = {" + ', '.join([str(l) for l in cont_bits_lengths]) + "};\n")

	program_source+=("	constant int comb_lengths["+str(len(comb_lengths))+"] = {" + ', '.join([str(l) for l in comb_lengths]) + "};\n")

	program_source+=("	constant int trajectory["+str(len(comb_lengths))+"] = {" + ', '.join([str(l) for l in t]) + "};\n")

	#n = nubmer of expected combinations in total
	#cont_bit = bit that is checked for


	expCombs = np.product([2**(i-1) for i in comb_lengths[a:b]])
	# print("expCombs", expCombs)
	program_source += 	"""
		kernel void stabiter(global long *a, global long *b){
			uint gid = get_global_id(0);
			long u = ((long)1)<<"""+str(qubits-1)+""";
			int l = """+str(int(expCombs))+""";
			long results[2]["""+str(int(expCombs))+"""];
			results["""+str(a%2)+"""][0] = a[gid];
			int r = 1;
			int n = 1;


			for (int bit = """+str(a)+"""; bit<"""+str(b)+"""; bit++){

				int length = comb_lengths[bit];

				int c = comb_lengths[bit];
				int cmm = comb_mask_mask[c];	

				int index = 0;

				for (int i=0; i<r;i++){
					long temp = results[(bit)%2][i];
					int bit_set = trajectory[bit];
					for(int j=0; j<cont_bits_lengths[bit]; j++){
						long bit_mask = u>>(cont_bits[bit][j]);
						bit_set ^= (temp&bit_mask)&&u;
					}
					for(int j=0; j<cmm; j++){
						results[(bit+1)%2][index] = temp;
			 		 	for (int k=0; k<comb_mask[c][bit_set][j]; k++){
			 		 		results[(bit+1)%2][index] = results[(bit+1)%2][index]^(u>>comb[bit][bit_set][j][k]);
			 		 	}
			 		 	index++;
					}			
				}
			 	r*=cmm;
			}
			for (int i=0; i<l; i++){
				b[gid*l+i] = results["""+str(b%2)+"""][i];
			}

		}

		"""
			# for (int i=0; i<128; i++){
			# 	b[gid+i] = results[0][i];
			# }
	bpcache[key] = (program_source, expCombs)
	return program_source, expCombs



def run_segment(traj, start, finish):
	source, expCombs = build_program(start, finish)
	# for line in source.split("\n"):print(line)

	program_source_bind = cl.Program(ctx, source)
	program = program_source_bind.build()


	a = np.array(traj).astype(np.uint64)
	b = np.empty([len(traj), int(expCombs)]).astype(np.uint64)

	a_buffer = cl.Buffer(ctx, flags=cl.mem_flags.READ_ONLY, size=a.nbytes)
	b_buffer = cl.Buffer(ctx, flags=cl.mem_flags.WRITE_ONLY, size=b.nbytes)

	queue = cl.CommandQueue(ctx)

	input_tuples = ((a, a_buffer),)
	output_tuples = ((b, b_buffer),)
	run_ocl_kernal(queue, program.stabiter, (len(a),), input_tuples, output_tuples)

	return b.flatten()
import random
bpcache2 = {}
def build_program2(traj):
	masks = []
	mods = []
	for stabilisers in {itertools.combinations(range(stabilisersX), r) for r in range(stabilisersX+1)}:
		for combo in stabilisers:
			mod = 1
			mask = np.zeros(qubits).view(np.uint64)
			combo_key = np.zeros(stabilisersX).view(np.uint64)
			for s in combo:
				if tx_bits[s] == 1:
					mod *= -1
				mask = np.bitwise_xor(mask, stabiliser_rawXnp[s])
			# print(combo, mod)
			masks.append(str(int(mask.dot(po2_array2))))
			mods.append(str(int(mod)))
	# {print(bkey, b2_dict[bkey]) for bkey in b2_dict}
	nkeys = len(masks)


	l = set()
	while True:
		eigenvalue = traj[random.randint(0,len(traj))]
		for mask in masks:
			k = int(eigenvalue)^int(mask)
			if k not in l:
				l.add(k)
		if len(l) == 2**(stabilisersX+1):
			break
		elif len(l) > 2**(stabilisersX+1):
			exception("something wrong")





	program_source = """
		#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
		constant int mlength = """ + str(nkeys) + """;
		constant int eigenmap["""+str(len(l))+"] = {" + ", ".join([str(e) for e in sorted(list(l))]) + """};
		constant int mods["""+str(nkeys)+"] = {" + ", ".join(mods) + """};
		constant int masks["""+str(nkeys)+"] = {" + ", ".join(masks) + """};"""


	program_source += 	"""
		kernel void xproc(global long *a, global long *b){
			uint gid = get_global_id(0);
			long eigenvalue = a[gid];
			int j = 0;

			for (int i = 0; i<""" + str(int(qubits)) + """; i++){
				if ((eigenvalue>>i)&1) {
					j++;
				}
			}

			for (int i = 0; i<"""+str(nkeys)+"""; i++){
				long new_eigenvalue = eigenvalue^masks[i];
				int memloc = 0;
				for (int k = 0; k<"""+str(len(l))+"""; k++){
					if (eigenmap[k] == new_eigenvalue){
						memloc = k;
						break;
					}
				}

				atom_add(&b[(memloc)*""" + str(qubits+1) +""" + j],mods[i]);	
			}

		}"""

	return program_source, l


# 	key = str(a)+"|"+str(b)
# 	if key in bpcache2:
# 		program_source, expCombs = bpcache2[key]
# 		return program_source, expCombs
# 	program_source = """

# 		constant int comb_mask_mask[4] = {0, 1, 2, 4};
# 		constant int comb_mask[4][2][4] = {
# 			{{-1, -1, -1, -1},{-1, -1, -1, -1}},
# 			{{0, -1, -1, -1}, {1, -1, -1, -1}},
# 			{{0, 2, -1, -1},{1, 1, -1, -1}},
# 			{{0, 2, 2, 2},{1, 1, 1, 3}}
# 		};
# 	"""




# 	max_combs = 0
# 	program_source+=("	constant int comb["+str(len(combinations))+"][2][4][3] = {\n") 
# 	first0 = True
# 	for i in range(len(combinations)):
# 		pad = np.zeros([2, 4, 3]).view(np.uint64)
# 		for j in [0,1]:
# 			for k in range(len(combinations[i][j])):
# 				for l in range(len(combinations[i][j][k])):
# 					pad[j][k][l] = combinations[i][j][k][l]
# 		first1 = True
# 		if first0:
# 			first0 = False
# 			program_source += "\n{"
# 		else:
# 			program_source += ",\n{"
# 		for i in pad:
# 			first2 = True
# 			if first1:
# 				first1 = False
# 				program_source += "\n{"
# 			else:
# 				program_source += ",\n{"
# 			for j in i:
# 				first3 = True
# 				if first2:
# 					first2 = False
# 					program_source += "\n	{"
# 				else:
# 					program_source += ",\n	{"
# 				program_source += ", ".join([str(k) for k in list(j)])
# 				program_source += "}"
# 			program_source += "}"


# 		program_source += "}"
# 	program_source += "};\n"




# 	cont_bits_lengths = []
# 	program_source += ("	constant int cont_bits["+str(len(cont_bits))+"][2] = {")
# 	first = True
# 	for bit in range(len(cont_bits)):
# 		print(cont_bits[bit])
# 		if first:
# 			first = False
# 			program_source += "\n"
# 		else:
# 			program_source += ",\n"
# 		pad = [-1, -1]
# 		for i in range(len(cont_bits[bit])):
# 			pad[i] = cont_bits[bit][i]


# 		program_source += "{" + ', '.join([str(l) for l in pad]) + "}"
# 		cont_bits_lengths.append(len(cont_bits[bit]))

# 	program_source += ("};\n")

# 	program_source+=("	constant int cont_bits_lengths["+str(len(cont_bits_lengths))+"] = {" + ', '.join([str(l) for l in cont_bits_lengths]) + "};\n")

# 	program_source+=("	constant int comb_lengths["+str(len(comb_lengths))+"] = {" + ', '.join([str(l) for l in comb_lengths]) + "};\n")

# 	program_source+=("	constant int trajectory["+str(len(comb_lengths))+"] = {" + ', '.join([str(l) for l in t]) + "};\n")

# 	#n = nubmer of expected combinations in total
# 	#cont_bit = bit that is checked for


# 	expCombs = np.product([2**(i-1) for i in comb_lengths[a:b]])
# 	program_source += 	"""
# 		kernel void stabiter(global long *a, global long *b){
# 			uint gid = get_global_id(0);
# 			long u = ((long)1)<<"""+str(qubits-1)+""";
# 			int l = """+str(int(expCombs))+""";
# 			long results[2]["""+str(int(expCombs))+"""];
# 			results["""+str(a%2)+"""][0] = a[gid];
# 			int r = 1;
# 			int n = 1;


# 			for (int bit = """+str(a)+"""; bit<"""+str(b)+"""; bit++){

# 				int length = comb_lengths[bit];

# 				int c = comb_lengths[bit];
# 				int cmm = comb_mask_mask[c];	

# 				int index = 0;

# 				for (int i=0; i<r;i++){
# 					long temp = results[(bit)%2][i];
# 					int bit_set = trajectory[bit];
# 					for(int j=0; j<cont_bits_lengths[bit]; j++){
# 						long bit_mask = u>>(cont_bits[bit][j]);
# 						bit_set ^= (temp&bit_mask)&&u;
# 					}
# 					for(int j=0; j<cmm; j++){
# 						results[(bit+1)%2][index] = temp;
# 			 		 	for (int k=0; k<comb_mask[c][bit_set][j]; k++){
# 			 		 		results[(bit+1)%2][index] = results[(bit+1)%2][index]^(u>>comb[bit][bit_set][j][k]);
# 			 		 	}
# 			 		 	index++;
# 					}			
# 				}
# 			 	r*=cmm;
# 			}
# 			for (int i=0; i<l; i++){
# 				b[gid*l+i] = results["""+str(b%2)+"""][i];
# 			}

# 		}

# 		"""
# 			# for (int i=0; i<128; i++){
# 			# 	b[gid+i] = results[0][i];
# 			# }
# 	bpcache2[key] = (program_source, expCombs)
# 	return program_source, expCombs
po2_array = 1 << np.arange(stabilisersX)[::-1]
po2_array2 = 1 << np.arange(qubits)[::-1]

def run_segment2(traj):
	source, l = build_program2(traj)
	# {print(line) for line in source.split("\n")}
	# c=0
	# l = {}
	# for eigenvalue in traj:
	# 	eigenvalue_bin = np.unpackbits(np.array([eigenvalue], dtype='>i8').view(np.uint8))[-qubits:]
	# 	j = np.count_nonzero(eigenvalue_bin)
	# 	for bkey in b2_dict:
	# 		mod = 1
	# 		#print(eigenvalue_bin, b2_dict[bkey], mod)
	# 		mask = np.zeros(qubits).view(np.uint64)
	# 		for bk in bkey:
	# 			if tx_bits[bk] == 1:
	# 				mod *= -1
	# 			mask = np.bitwise_xor(mask,stabiliser_rawXnp[bk])
	# 		#print(mask)
	# 		k = np.bitwise_xor(eigenvalue_bin,mask).dot(po2_array2)
	# 		if k not in l:
	# 			l[k] = np.zeros(qubits+1)
	# 		#print(k,j,mod)
	# 		l[k][j] += mod

			
			
	# 	c += 1
	# 	#print("done",eigenvalue, c,"/"+str(len(traj)))
	# 	break
	# #print(l)
	# return l

	program_source_bind = cl.Program(ctx, source)
	program = program_source_bind.build()


	a = np.array(traj).astype(np.int64)
	b = np.empty([len(l), int(qubits+1)]).astype(np.int16)

	a_buffer = cl.Buffer(ctx, flags=cl.mem_flags.READ_ONLY, size=a.nbytes)
	b_buffer = cl.Buffer(ctx, flags=cl.mem_flags.WRITE_ONLY, size=b.nbytes)

	queue = cl.CommandQueue(ctx)

	input_tuples = ((a, a_buffer),)
	output_tuples = ((b, b_buffer),)
	run_ocl_kernal(queue, program.xproc, (len(a),), input_tuples, output_tuples)

	return list(sorted(l)), b

tt = time.time()
logic = np.zeros(2*qubits+2, np.uint64)
trajectories = run_segment(trajectories, 0, int(num_stabilisers/3))
trajectories = run_segment(trajectories, int(num_stabilisers/3), int(num_stabilisers*2/3))
trajectories = run_segment(trajectories, int(num_stabilisers*2/3), num_stabilisers)
keys, zzz = run_segment2(trajectories)
z = 0
# for i in range(len(keys)):
# 	print(keys[i], zzz[i])
# trajectories = run_segment(trajectories, 13, 18)
# print(sorted(trajectories))
print(time.time()-tt)
#cpu
#d2 0.09999
#d3 0.113
#d4 2.37
#d5 645.11

# clen = len(trajectories)
# n = 10000
# start = time.time()
# final_traj = np.zeros(2*qubits+2+210000).astype(np.uint64)
# for tr in range(0, clen, n):
# 	prelim = run_segment(trajectories[tr:tr+n], 18, 25)
# 	final_traj = run_segment2(prelim, 25, 30, final_traj)
# 	#print(list(final_traj[210000:qubits+210001]))
# 	p = round(tr/clen, 8)

# 	try:
# 		print(p, time.time()-start, (time.time()-start)*((1-p)/p))
# 	except:
# 		pass


# #10 21000 1074150
# #100 2100 10627
# #1000 2.23 1115.07
# #10000 0.58 300
# #100000 0.079 53.60
# #1000000 0.047


# print(list(final_traj[210000:qubits+210001]))
# print(list(final_traj[qubits+210001:]))

#traj
#traj as bits
#j = count 1s
#ktemp = eigenvalue
#loop through all combinations of X, e.g. 00, 01, 10, 11
	#mod = 1
	#if bit == 1 and x target==1, flip mod (e.g. targetX 11, combo 01, mod 1 -> -1)
	#flip any bit involved with the stab  (e.g. k_temp = 0000, stabiliser 012, k_temp = 1110)
	#logic[k(po2)][j]+=mod (where logic[k(po2)] starts as an array of 0s qubits+1 long)



# IIXXX
# XXXII
# target 01

# traj:
# ABCDE
# j = countones(ABCDE)

# mod = 1

# loop 00, 01, 10, 11
# mod 1, -1, 1, -1
# ABCDE, abcDE, ABcde, ABcDE

# l[ABCDE] += 1
# l[abcDE] -= 1
# l[ABcde] += 1
# l[ABcDE] -= 1


# eigenvalue = EFGHJ
# target = AB
# loop = CD
# mod = AB and CD bitwise, count bits
# mask = A(stab)if1 xor B(stab)if1

# l[eigenvalue^mask][j] += mod