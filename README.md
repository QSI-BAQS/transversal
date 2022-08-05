# Transveral Injection - Encoding Non-Clifford States Using The Surface Code 

Transversal injection is a novel method of preparing non-clifford ancillary states. This is important for realising non-Clifford gates for universal quantum computation. This code repo contains:
a) A rotated surface code full state simulation
b) A python wrapper to benchmark the surface code

## Requirements
gcc
python 3.9 or later (earlier version can be used if lru_caching is disabled in `transversal.py` and `stabiter.py`)

## Setup and Usage
1) Compiler the simulated using `gcc main.c Simulator/sim.c Simulator/norm.c -lm -o sim`

2) Configure simulation
- See next section

3) Run the experiment either using
- command line execution `python3.9 transversal.py` 
*or* 
- module import and function call `import transversal` followed by `run_batch(distance, ti_flag, angle_theta, angle_phi, err, fast_mode, fn, reps)`


## Configuration
To configure the experiment, the following values can be changed to suit. 

```Python
	distance # integer: surface code distance
	ti_flag # boolean: enables transversal injection, if disabled this is essentially just the surface code
	angle_theta # float: theta in radians 
	angle_phi # float: phi in radians 
	err # float: physical error rate
	fast_mode # boolean: enabling fast_mode will only run full simulations when an error occurs. This is significantly faster but not suitable for post-select experiments. 
	fn # string: file name for storing results
	reps # integer: if fast_mode disabled number of repetitions - else, number of experiments with errors
```

## Results
Results of the simulation will be displayed during execution. A dictionary containing results on a per-trajectory basis will be returned. A backup will be saved in /results/ every 1000 samples.