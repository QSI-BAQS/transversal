import numpy as np
from functools import cache, lru_cache
import rotated_layout
import itertools
import gc
global loopsdistance

def commute(stabs, target):
	commutes = []
	global loops
	t0 = np.zeros(stabilisersZ)
	cc = 0

	for c in range(1, len(stabs)+1):
		for comb in itertools.combinations(stabs,c):
			stabcomb = d[np.array(comb)]
			xor = np.bitwise_xor.reduce(stabcomb)
			if (xor==target).all():
				commutes.append((comb, xor))
			cc += 1
			loops += 1


	if (commutes == []):
		return False
	return commutes


@lru_cache(maxsize=1000000)
def findones(args):
	scache, tcache = args
	if scache:
		availablestabs = set([int(c) for c in scache.split('|')])
	else:
		return False

	target = np.unpackbits(np.array([tcache], dtype='>i8').view(np.uint8))[-stabilisersZ:]

	results = []
	if (target == t0).all():

		tmasked = target
		stabs = set(availablestabs)
		nresults = []
		for stab in availablestabs:
			ntarget = np.array(d[stab])
			nstabs = set(availablestabs)
			nstabs.remove(stab)
			nresults = findones(cacheopt(nstabs, ntarget))
			if nresults:
				for r in nresults:
					results.append([stab] + list(r))

	else:
		tmasked = np.ma.array(target, mask = np.bitwise_xor(target, np.ones(len(target), dtype=int)))
		stabs = set()
		for i in range(len(tmasked)):
			if tmasked[i]:
				for s in stabdirectory[i]:
					if s in availablestabs:
						stabs.add(s)
		if stabs == set():
			return False


		candidatecombs = commute(stabs, tmasked)
		if candidatecombs:
			for comb in candidatecombs:
				newstabs = set(availablestabs)
				newstabs.difference_update(stabs)
				if (comb[1] == target).all():
					results.append(comb[0])
				if newstabs != set():
					new_target = np.bitwise_xor(comb[1], target)
					nresults = findones(cacheopt(newstabs, new_target))
					if nresults:
						for r in nresults:
							results.append(list(comb[0]) + list(r))
	return results

def cacheopt(s, t):
	return ('|'.join(str(i) for i in s), t.dot(po2_array))






def algorithm_2(target, stabdirectory, qubits): # python implementation of algorithm 2, opencl implementaiton can go to d8 and higher but in experimental state
	iworkspace = []
	farray = np.array([]).astype(int)
	fworkspace = []
	usedstabs = set()
	combs_record = ""
	for i in range(len(target)):
		ll = 0
		if len(fworkspace)>1000000:
			farray = np.append(farray, np.array(fworkspace).astype(int))
			while fworkspace:
				g = fworkspace.pop()
			del fworkspace   # desperate attempts at managing memory in python - eww
			fworkspace = []
			gc.collect()
			gc.collect(0)
			gc.collect(1)
			gc.collect(2)  

		if len(iworkspace) == 0: #First run
			for c in range(target[i], len(stabdirectory[i])+1, 2):
				for comb in itertools.combinations(stabdirectory[i], c):
					traj = 0
					for q in comb:
						traj |= 1<<(qubits-q-1)
					fworkspace.append(traj)

		else:
			for traj in iworkspace:
				ti = target[i]
				unused = set(stabdirectory[i])
				for q in stabdirectory[i]:
					if q in usedstabs:
						ti ^= (traj&(1<<(qubits-q-1))>0)
						unused.remove(q)

				ll = len(unused)
				for c in range(ti, len(unused)+1, 2):
					for comb in itertools.combinations(unused, c):
						ntraj = int(traj)
						for q in comb:
							ntraj |= 1<<(qubits-q-1)
						fworkspace.append(ntraj)
		

		combs_record += str(ll)
		del iworkspace

		iworkspace = np.append(farray, np.array(fworkspace).astype(int))
		while fworkspace:
			g = fworkspace.pop()
		del fworkspace

		del farray

		fworkspace = []
		farray = np.array([]).astype(int)

		gc.collect()
		gc.collect(0)
		gc.collect(1)
		gc.collect(2)
		for q in stabdirectory[i]:
			usedstabs.add(q)
	return iworkspace

def solve(distance, stab_int):
	qubits = distance**2

	stabiliser_rawX, stabiliser_rawZ = rotated_layout.get_stabs(distance, True)
	stabilisersX = len(stabiliser_rawX)
	stabilisersZ = len(stabiliser_rawZ)


	stabiliser_rawX = list(reversed(stabiliser_rawX))
	stabiliser_rawZ = list(reversed(stabiliser_rawZ))

	loops = 0
	po2_array = 1 << np.arange(stabilisersZ)[::-1] #Used to make the conversion of bin array to int fast
	po2_array2 = 1 << np.arange(qubits)[::-1] #Used to make the conversion of bin array to int fast



	stabiliser = np.array(stabiliser_rawZ).astype(int)

	stabdistil = []
	stabdistil2 = []
	for i in range(qubits):
		stabdistil.append(np.count_nonzero(stabiliser==i, axis=1))
		stabdistil2.append(np.count_nonzero(stabiliser==i, axis=1).dot(po2_array))

	t0 = np.zeros(stabilisersZ).astype(int)

	d = np.array(stabdistil)

	tx = np.unpackbits(np.array([stab_int], dtype='>i8').view(np.uint8))[-stabilisersX-stabilisersZ:-stabilisersZ]
	tz = np.unpackbits(np.array([stab_int], dtype='>i8').view(np.uint8))[-stabilisersZ:]

	stabs = set(range(qubits))

	path0 = [[],tz]
	stabdirectory = [set() for i in range(stabilisersZ)]
	for i in range(len(stabdistil)):
		for qubit in range(len(stabdistil[i])):
			if stabdistil[i][qubit] == 1:
				stabdirectory[qubit].add(i)




	combs2 = algorithm_2(tz, stabdirectory, qubits)

	done = set()
	logic = {}



	for traj in combs2: # This loop looks at the Z-trajectories that have been considered valid for the stabiliser measurements and now projects into the X stabilisers before sorting it into the final data structure
		traj_array = np.unpackbits(np.array([traj], dtype='>i8').view(np.uint8))[-qubits:]
		j = np.count_nonzero(traj_array)
		for stabilisers in {itertools.combinations(range(stabilisersX), r) for r in range(stabilisersX+1)}:
			for combo in stabilisers:
				mod = 1
				k_temp = np.array(traj_array)
				for s in combo:
					if s != ():
						if tx[s] == 1:
							mod *= -1
						for bit in stabiliser_rawX[s]:
							if bit != -1:
								k_temp[bit] ^= 1

				k = k_temp

				try:
					logic[k.dot(po2_array2)][j] += mod #
				except:
					logic[k.dot(po2_array2)] = np.zeros(qubits+1, int)
					logic[k.dot(po2_array2)][j] += mod
	eqs = []
	for key in logic:
		skey = list(logic[key])
		if skey in eqs:
			logic[key] = eqs.index(skey)
		else:
			index = len(eqs)
			eqs.append(skey)
			logic[key] = index


	
	return logic, eqs



### How to use this module:
# stab_int = int('00110100',2) # convert the binary string representation of stabiliser measurements to an integer, i.e. in the format {X_1..n, Z_1...n}
# logic = solve(3, stab_int) # call solve(distance, stabiliser_integer)
