#!/bin/env python3
'''
Expand all files from the case with the specified ID number in ./history, overwriting.

Usage: in run directory, `rest.py <id number e.g. 5>`
'''
import os
import sys
import subprocess
import glob

try:
    num = sys.argv[1]
except:
    print('rest <archive number>')
    sys.exit()
    
if 'SOLPSTOP' not in os.environ.keys():
    print('Error: SOLPS not initialized. Exiting.')
    sys.exit(0)
    
files = glob.glob(f'history/{num}/20*.tar.gz')
if len(files) != 1:
    print('Invalid archive number')
    sys.exit()
filename = files[0]

subprocess.call(['tar', '-xf', filename])

if os.path.isfile('b2fplasmf'):
    print('Removing b2fplasmf')
    os.remove('b2fplasmf')
    if os.path.isfile('b2uf.prt'):
        os.remove('b2uf.prt')
subprocess.call(['myb2uf'])
