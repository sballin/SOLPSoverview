#!/bin/env python3
'''
Delete run ID files, including start and finish.
'''
import sys
import os
import shutil

try:
    if len(sys.argv) == 2:
        startid = int(sys.argv[1])
        endid = int(sys.argv[1])
    elif len(sys.argv) == 3:
        startid = int(sys.argv[1])
        endid = int(sys.argv[2])
except:
    print('del.py <start id> <end id (optional)>')
    sys.exit()
    
to_remove = []

for i in range(startid, endid+1):
    to_remove.append(f'history/{i}')

print('Will remove the following files:')
for f in to_remove:
    print(f)
    
answer = ''
while answer not in ["y", "n"]:
    answer = input("Confirm? [y/n] ").lower()
if answer == "y":
    for file in to_remove:
        try:
            shutil.rmtree(file)
        except Exception as e:
            print(file, ':', e)

