#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 11:47:29 2017

@author: yboetz
"""

import os, sys, argparse
import numpy as np
import subprocess as sp


parser = argparse.ArgumentParser(description='Scaling of hydro code')
parser.add_argument('-g', nargs = '+', dest="g", help = 'Gridsize', type = int)
parser.add_argument('-s',  action='store', dest="s", help = 'Steps', type = int)
parser.add_argument('-n', nargs = '+', dest="n", help = 'Number of CPUs per socket', type = int)
parser.add_argument('--omp', action="store_true", help = 'Use OMP')
parser.add_argument('--hwthread', action="store_true", help = 'Use MPI hwthreads')
options = parser.parse_args(sys.argv[1:])

nx, ny = options.g if options.g else (300,100)
steps = options.s if options.s else 1000
n = options.n if options.n else 1
threads = 2 if options.omp else 1
bind = 'hwthread' if options.hwthread else 'core'

inputText = """This namelist contains various input parameters for HYDRO runs

&RUN
nstepmax={steps}
tend=5000.0
noutput=10000000
on_output=.true.
/

&MESH
nx={nx}
ny={ny}
dx=0.01
boundary_left=1
boundary_right=1
boundary_down=1
boundary_up=1
idimbloc=21
jdimbloc=21
/

&HYDRO
courant_factor=0.8
niter_riemann=10
/
"""

rand = np.random.randint(10**12)
inFile = "input_" + str(rand)
outFile = "output_" + str(rand)

time = []
for cpu in n:
    myenv = os.environ.copy()
    myenv["OMP_NUM_THREADS"] = str(threads)

    keys = {'steps': steps, 'nx': nx, 'ny': ny, 'cpu': cpu, 'ifile': inFile, 'ofile': outFile, 'bind': bind}
    with open(inFile,'w') as file, open(outFile,'w') as _:
        file.write(inputText.format(**keys))
    
    cmd = 'mpirun -np {cpu} --map-by {bind} --bind-to {bind} hydro_mpiomp -i {ifile} > {ofile}'.format(**keys)
    try:
        sp.call(cmd, shell=True, env = myenv)
    except Exception as e:
        print(e)
    
    with open(outFile,'r') as file:
        lines = file.readlines()
        line = lines[-1]
        print(lines[-2],end='')
        begin, end = line.find('('), line.find(')')
        try:
            time.append(float(line[begin+1:end]))
        except Exception as e:
            print(e)

print(time)

try:
    os.remove(inFile)
    os.remove(outFile)
except Exception as e:
    print(e)