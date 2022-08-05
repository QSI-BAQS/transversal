# generates the stabilisers for distance d rotated surface code
# padded indicates whether or not the result should be padded with -1, only relevant for stabiter function
def get_stabs(d, padded):
	grid = [[0 for i in range(2*d+1)] for j in range(2*d+1)]
	c = 0
	for j in range(1, d*2+1, 2):
		for i in range(1, d*2+1, 2):
			grid[j][i] = ("L", c)
			c+=1

	c = 0
	for j in range(0, d*2+1, 2):
		for i in range(2, d*2-1, 2):
			if (i%4==2 and j%4==0) or (i%4==0 and j%4==2):
				grid[j][i] = ("X", c)
				c+=1

	c = 0
	for j in range(2, d*2-1, 2):
		for i in range(0, d*2+1, 2):
			if (i%4==0 and j%4==0) or (i%4==2 and j%4==2):
				grid[j][i] = ("Z", c)
				c+=1

	ancillas_X = []
	ancillas_Z = []


	for j in range(2*d+1):
		for i in range(2*d+1):
			if grid[j][i] != 0:
				if grid[j][i][0] == "X":
					ancillas_X.append([])
					for a in [1,-1]:
						for b in [1,-1]:
							try:
								ancillas_X[-1].append(grid[j+a][i+b][1])
							except:
								pass
					if padded:
						while len(ancillas_X[-1]) != 4:
							ancillas_X[-1].append(-1)

				if grid[j][i][0] == "Z":
					ancillas_Z.append([])
					for a in [1,-1]:
						for b in [1,-1]:
							try:
								ancillas_Z[-1].append(grid[j+a][i+b][1])
							except:
								pass
					if padded:
						while len(ancillas_Z[-1]) != 4:
							ancillas_Z[-1].append(-1)

	return ancillas_X, ancillas_Z

