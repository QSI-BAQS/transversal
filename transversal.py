import math
import subprocess
import cmath
import numpy as np
import os
import time
import stabiter
from functools import cache, lru_cache



# Searches for 'solverhelper.csv' - contains a precomputed dictionary produced from algorithm 1 
try:
	f = open("solverhelper.csv", "r")
	dict_found = 4 # set to the distance contained in the file
	fsize = os.stat('solverhelper.csv').st_size
except:
	dict_found = False
	fsize = -1
	f = -1


@lru_cache(maxsize=1000000)
def geteq(d,s):
	if (not dict_found) or d!=dict_found: 
		return stabiter.solve(d, s)
	else:
		if s>3000: # kinda smart kinda dumb algorithm for searching through a long file quickly
			f.seek(int(fsize*s/(2**15)*0.8))
			i = -2
		else:
			f.seek(0)
			i = 0
		for line in f:
			if i == -1:
				# print(line, line.split(' ')[0])
				i = int(line.split(' ')[0])
			elif i == s:
				eqn = eval('('+line.split('(')[-1].strip())
				return eqn
			i += 1



def polar(x, y):
	return cmath.polar(x), cmath.polar(y)


def bloch(x, y):
	if not x and y:
		return (0, complex(1, 0))
	if not y and x:
		return (1, complex(0, 0))
	if not y and not x:
		return (0, complex(0, 0))
	a, b = polar(x, y)
	b = cmath.rect(b[0],b[1] - a[1])
	a = a[0]
	mag = np.sqrt(a**2 + b.real**2 + b.imag**2)
	if mag == 0:
		mag = 1
	return a/mag, b/mag

# Computes expected final state based on initial transversal rotation of th radians in Y direction and the given analytical equation for that trajectory
def guess(th, ph, d, xz, ti):  
	if ti:
		a, b = initial(th, ph)
	else:
		a, b = (1, complex(0, 0))

	eqs = geteq(d,xz)

	eq0 = np.array(eqs[1][0])
	terms = len(eqs[1][0])
	results = []

	for eq in eqs[1]:
		z = sum([eq[i]*pow(a, terms-i)*pow(b, i) for i in range(len(eq))])
		results.append(z)

	svdict = {}
	for stateint in eqs[0]:
		svdict[stateint] = results[eqs[0][stateint]]
	return svdict

# Remove the ancilla from the raw simualtion results
def removeAncilla(stateout, distance):
	states = {}
	stateskey = []
	states_actual = []
	svdict = {}
	for line in stateout:
		if line.strip():
			spl = line.split(b":")
			stateint = int(bin(int(spl[0]))[2:-1].zfill(distance**2),2) # bitshift to account for ancilla, >>1 was being annoying
			if abs(float(spl[1].strip()))<0.00000001:
				x = 0
			else:
				x = float(spl[1].strip())

			if abs(float(spl[2].strip()))<0.00000001:
				y = 0
			else:
				y = float(spl[2].strip())

			z = complex(x, y)
			svdict[stateint] = z		
	return svdict

# Performs a tensor product of the simulated and reference state to assess fidelity
def compareStates(g, a, distance):
	sv_g = np.array(np.zeros(2**(distance**2)), np.dtype(np.complex64))
	sv_a = np.array(np.zeros(2**(distance**2)), np.dtype(np.complex64))
	if distance == 3:
		for i in g:
			mod = 1
			if bin(i)[2:].zfill(9)[:3].count('1')%2 == 1:
				mod = -1
			sv_g[i] = (g[i])*mod

		for i in a:
			mod = 1
			if bin(i)[2:].zfill(9)[:3].count('1')%2 == 1:
				mod = -1
			sv_a[i] = (a[i])
	else:
		for i in g:
			sv_g[i] = (g[i])

		for i in a:
			sv_a[i] = (a[i])

	norm_g = np.linalg.norm(sv_g)
	if norm_g: sv_g = sv_g/norm_g


	norm_a = np.linalg.norm(sv_a)
	if norm_g: sv_a = sv_a/norm_a

	result = abs(np.tensordot(sv_a, sv_g,axes=1).item())
	return result

# Returns the initial state of a qubit rotated by th radians about the Y axis. This can be replaced to start at an arbitrary alpha/beta although the simulator file will need to be changed to reflect this
def initial(th, ph):
	return math.cos(th/2), complex(math.cos(ph), math.sin(ph))*math.sin(th/2)


def run_sim(err: float, dist: int, angle_theta: float, angle_phi: float, ti: bool, round_len:int, fast_mode):
	empty = True
	runs = 0
	while True: 
		x = subprocess.run(["./sim", str(dist), str(err), ("1" if ti else "0"), str(angle_theta), str(angle_phi), ("1" if fast_mode else "0")], stdout=subprocess.PIPE)
		if x.returncode==0: #run until success, simulation will fail upon error detection
			break
		else:	
			runs += int(x.stdout.split(b'\n')[2].split(b' ')[-1])


	out = x.stdout.split(b'\n')
	errors = {'X':0, 'Z':0, 'XZ':0}
	statesect = False
	stateout = []
	rawout = []
	staberr = []
	i = 0
	t = 0
	for line in out:
		if b"Runs:" in line:
			runs += int(line.split(b' ')[-1])
		if b"abort" in line:
			empty = True
			return "aborted", empty, stateout, rawout, runs
		if statesect:
			if b'state' not in line:
				stateout.append(line)
		else:
			if b'state vector:' in line:
				statesect = True
			if b'bit' in line:
				t = math.floor(i/round_len)
				i+=1
				empty = False
				spl = line.split(b' ')
				rawout.append(int(spl[-1]))


	return errors, empty, stateout, rawout, runs



def run_batch(distance, ti, angle_theta, angle_phi, err, fm, fn, reps):
	trajresults = {}
	runtime = time.time()
	global fast_mode


	stabs = stabiter.rotated_layout.get_stabs(distance,False)
	round_len = len(stabs[0])+len(stabs[1])
	repsc = reps
	right = 0
	wrong = 0
	precomp = {}
	while (reps > 0):

		# reports progress
		if (distance==2 and (reps%100 == 0)) or (distance>2 and reps%10 == 0):
			pcd = ((repsc-reps)/repsc) + 10**-10
			runningfor = time.time()-runtime
			pctg = 1 - pcd
			if reps<repsc:
				print(str(round(pcd*100,2)) + "% Time elapsed:", round(runningfor/60,2), "minutes, time to go:", round(runningfor/pcd*pctg/60,2))

		# record results
		if (reps%1000 == 0):
			try:
				with open("results/"+fn, 'w') as f_b:
					f_b.write(str(trajresults))
					print("backed up")
			except:
				os.mkdir("results")

		reps -= 1

		errors, empty, stateout, raw, runs = run_sim(err, distance, angle_theta, angle_phi, ti, round_len, fm)
		i = 0 

		states = 0
		for i in stateout:
			if abs(float(stateout[1].split(b': ')[0]))>0.00000001:
				states += 1


		if not empty:
			raw_cat = ''.join([str(i) for i in raw[:round_len]])
			raw_fat = ''.join([str(i) for i in raw])
			traj = int(raw_cat,2)


			right += runs				


			if traj not in precomp:
				g = guess(angle_theta, angle_phi, distance, traj, ti)
				precomp[traj] = g
			else:
				g = precomp[traj]

			

			a = removeAncilla(stateout, distance)

			
			fidelity = compareStates(g, a, distance)

			if traj not in trajresults:
				trajresults[traj] = {"runs":runs, "right":0, "wrong":0}
			else:
				trajresults[traj]["runs"] += runs

			if fidelity > (1-err-0.00001):
				right += 1
				result = "right"
			else:
				wrong += 1
				result = "wrong"
				# break

			trajresults[traj][result] += 1

			if (distance==2 and (reps%100 == 0)) or (distance>2 and reps%10 == 0):
				print("distance", distance, "err_Physical:", err, "right", int(right), "wrong", int(wrong), "err_Logical", wrong/(right+wrong))


	print("Total runtime:", time.time()-runtime)
	return trajresults





if __name__ == '__main__':
	# Code distance
	distance = 3
	# # Transversal injection ON 
	ti_flag = True 
	# # Transversal rotations
	angle_theta = math.pi*45/180 
	angle_phi = math.pi*45/180 
	# # Physical error rate
	err = 0.0001
	# # Skips runs without errors, not suitable for post-selection data
	fast_mode = True 
	# # Results filename
	fn = str(distance) + "_" + str(err)[2:] + "_" + str(time.time())
	# # Number of experiments (Number of errored experiments if fast_mode is on)
	reps = 5000 

	run_batch(distance, ti_flag, angle_theta, angle_phi, err, fast_mode, fn, reps)
