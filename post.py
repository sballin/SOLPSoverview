#!/bin/env python3
'''
After a run is over, clean up prt files, make plots and point to case .tar.gz, saving files in ./history under the user-specified ID number.

Usage: in run directory do e.g. `post.py 5` for the run with ID 5 
'''
import sys
import shutil
import subprocess
import glob
import os

def get_latest(file_searchstring):
    files = glob.glob(file_searchstring)
    return max(files, key=os.path.getctime)

with open('current_run.txt', 'r') as f:
    current_run = f.read()
    num = int(current_run.split()[0].split('id')[1])
    print(f'Doing post id{num}')
    slurmid = int(current_run.split()[1])

run_failed = os.path.exists('b2mn.exe.dir') and not os.path.isfile('b2fstate')

if run_failed:
    print('Run failed')
    env = dict(os.environ)
    files = glob.glob('b2mn.exe.dir/plasmastate*')
    for f in files:
        env['PLASMASTATE'] = f
        suffix = f.split('plasmastate.')[-1]
        subprocess.call(['b2run', 'b2co'], env=env)

        print(f'Copying b2fstate to history/b2fstate{num}_{suffix}')
        shutil.copyfile('b2fstate', f'history/b2fstate{num}_{suffix}')

        print('Making plots')
        subprocess.call(['plotall'])
        print(f'Copying plots.pdf to history/plots{num}_{suffix}.pdf')
        shutil.copyfile('plots.pdf', f'history/plots{num}_{suffix}.pdf')
else:
    try:
        os.remove('b2mn.prt')
    except:
        print('No b2mn.prt to remove')
        
    try:
        os.mkdir(f'history/{num}')
    except Exception as e:
        pass
        
    for file in [f'SLURM-{slurmid}.out', f'SLURM-{slurmid}.err']:
        try:
            shutil.move(file, f'history/{num}/{file}')
        except Exception as e:
            print(e)

    print('Making plots')
    subprocess.call(['plotall'])
    subprocess.call(['chmod', '644', 'plots.pdf']) # started being 600 17.2.25
    print(f'Copying plots.pdf to history/{num}/plots{num}.pdf')
    shutil.copyfile('plots.pdf', f'history/{num}/plots{num}.pdf')
    
    print('Doing quicksave')
    subprocess.call(['quicksave'])
    print(f'Moving quicksave.h5 to history/{num}/quicksave{num}.h5')
    shutil.move('quicksave.h5', f'history/{num}/quicksave{num}.h5')

    print('Saving run settings')
    subprocess.call(['save_run_settings', f'history/{num}/run_settings{num}.txt'])
    
    print('Making table')
    subprocess.call(['make_table'])
    
    latest_file = get_latest('20*.tar.gz')
    shutil.move(latest_file, f'history/{num}/{latest_file}')
    latest_toc = get_latest('20*.tar.toc')
    shutil.move(latest_toc, f'history/{num}/{latest_toc}')

    print('Uploading PDF and table')
    path = os.getcwd().split('/')
    ret1 = subprocess.call(['ssh', '-t', 'mfe', 'mkdir', '-p', f'~/public_html/solps/{path[-2]}/{path[-1]}'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    plotfile = f'history/{num}/plots{num}.pdf'
    ret2 = subprocess.call(['scp', plotfile, f'mfe:~/public_html/solps/{path[-2]}/{path[-1]}'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if ret1 != 0 or ret2 != 0:
        print('Error uploading PDF')
    ret3 = subprocess.call(['scp', 'history/index.html', f'mfe:~/public_html/solps/{path[-2]}/{path[-1]}'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if ret3 != 0:
        print('Error uploading table')
    print(f'https://www1.psfc.mit.edu/~sballinger/solps/{path[-2]}/{path[-1]}')
