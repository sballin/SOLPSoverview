#!/bin/env python3
'''
Do multiple runs with the same settings, incrementing the ID automatically from the largest one in ./history, saving plots PDF, and updating b2fstati after each completion.

Usage: in run directory, `run.py <number of runs e.g. 10>`
'''
import sys
import signal
import shutil
import random
import subprocess
import glob
import os
import h5py
import time
from datetime import timedelta
import argparse
from collections import defaultdict
from rich.progress import Progress, TextColumn, BarColumn, MofNCompleteColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn, Text, ProgressColumn

def sigint_handler(sig, frame):
    print(f'scancel {jobid}')
    subprocess.call(['scancel', jobid])
    sys.exit(0)
    
def sigquit_handler(sig, frame):
    global userquit, refresh_seconds
    print('touch b2mn.exe.dir/.quit')
    subprocess.call(['touch','b2mn.exe.dir/.quit'])
    userquit = True
    refresh_seconds = 10

def get_latest_save():
    files = glob.glob('history/[0-9]*')
    if len(files) == 0:
        return 0
    numbers = []
    for f in files:
        try:
            numbers.append(int(f.split('history/')[1]))
        except:
            pass
    return max(numbers)

def job_status(jobid):
    '''
    Output examples:
    '   PENDING ': waiting to run
    '   RUNNING \n   RUNNING \n   RUNNING ': SOLPS running
    ' COMPLETED \n COMPLETED \n COMPLETED ': SOLPS finished cleanly before job time limit
    '   TIMEOUT \n CANCELLED \n COMPLETED ': SOLPS was still running and hit job time limit
    Error output examples:
    'sacct: fatal: Bad job/step specified: --format=State'
    ''
    '''
    try:
        # --allocations avoids showing 3 rows per job
        out = subprocess.getoutput(f'sacct --job {jobid} --format=State --noheader --allocations')
        if not 'PENDING' in out and not 'RUNNING' in out and not 'COMPLETED' in out:
            print('Unusual job status:', out)
        return out.strip()
    except Exception as e:
        print('job_status exception:', e)
        return 'error'

def is_in_queue(jobid):
    global notified_squeue_err
    try:
        out = subprocess.getoutput('squeue -u $USER')
        if 'squeue: error' in out:
            # notify and keep running if squeue: error: Invalid user: $USER
            if not notified_squeue_err:
                try:
                    send_email(f'squeue error', out)
                    notified_squeue_err = True
                except Exception as e:
                    print('send_email about squeue error failed:', e)
            print(out)
            return True
        for line in out.split('\n'):
            if jobid in line and 'node' in line:
                return True
    except Exception as e:
        print('is_in_queue exception:', e)
    return False

def getIterations():
    if os.path.exists("b2mn.exe.dir/b2ftrace"):
        filename = "b2mn.exe.dir/b2ftrace"
    elif os.path.exists("cleanup/b2ftrace"):
        filename = "cleanup/b2ftrace"
    elif os.path.exists("b2ftrace"):
        filename = "b2ftrace"
    else:
        return 0
    with open(filename) as f:
        lines = f.readlines()

    mydatalist = []
    counter = 0
    read_data = False
    for line in lines:
        if "data" in line:
            counter += 1
            read_data = True
            continue
        if read_data:
            line = line.split()
            part_list = [float(i) for i in line]
            mydatalist.extend(part_list)
        
    return counter

def parseDat(filename):
    out = {}
    with open(filename, 'r') as f:
        txt = f.read()
    lines = txt.split('\n')
    for l in lines:
        try:
            if l[0] == "'":
                l = l.split('#')[0]
                key, val = l.split() 
                key = key.replace("'", "")
                val = float(val.replace("'", ""))
                out[key] = val
        except Exception as e:
            pass
    return out

def send_email(subject, body):
    recipient = 'sballin@mit.edu' #subprocess.getoutput('get_git_email') # in solps-iter/scripts
    # Create the email content
    message = f"""Subject: {subject}
To: {recipient}

{body}
"""
    # Write the message to a temporary file
    with open('email.txt', 'w') as file:
        file.write(message)
    # Send the email using sendmail
    subprocess.run(['/usr/sbin/sendmail', recipient], input=message.encode())
    os.remove('email.txt')

class TimeRemainingColumnCustom(ProgressColumn):
    """Renders estimated time remaining.

    Args:
        compact (bool, optional): Render MM:SS when time remaining is less than an hour. Defaults to False.
        elapsed_when_finished (bool, optional): Render time elapsed when the task is finished. Defaults to False.
    """

    def __init__(
        self,
        compact = False,
        elapsed_when_finished = False,
        table_column = None,
    ):
        self.compact = compact
        self.elapsed_when_finished = elapsed_when_finished
        super().__init__(table_column=table_column)

    def render(self, task: "Task") -> Text:
        """Show time remaining."""
        if self.elapsed_when_finished and task.finished:
            task_time = task.finished_time
            style = "progress.elapsed"
        else:
            task_time = task.time_remaining
            style = "progress.remaining"

        if iterations_completed < 0:
            return Text("-:--", style=style)
        
        minutes, seconds = divmod(int(seconds_remaining), 60)
        hours, minutes = divmod(minutes, 60)

        return Text(f"-{hours:d}:{minutes:02d}", style=style)

def fail_info():
    out = ''

    # Show job info to identify bad nodes
    try:
        out += subprocess.getoutput(f'sacct -j {jobid} --format=JobID,JobName,Elapsed,State,NodeList') + '\n\n'
    except Exception as e:
        print('fail_info: get jobinfo fail:', e)

    # Look for state which indicates run failed for b2/eirene numerics reasons
    if os.path.exists('b2mn.exe.dir') and not os.path.isfile('b2fstate'):
        out += 'Run failed: b2mn.exe.dir exists and b2fstate does not.\n\n'

    # Find and read SLURM output files
    if os.path.exists(f'history/{sim}/SLURM-{jobid}.out'):
        outstub = f'history/{sim}/SLURM-{jobid}'
    elif os.path.exists(f'SLURM-{jobid}.out'):
        outstub = f'SLURM-{jobid}'
    else:
        out += f'SLURM-{jobid}.out not found'
        return out
    with open(outstub+'.out', 'r') as f:
        outtxt = f.read()
    with open(outstub+'.err', 'r') as f:
        errtxt = f.read()

    # Display lines containing string 'err' and how often they repeated
    out += f'All lines of {outstub}.err/out containing err:\n\n'
    alltxt = errtxt+'\n'+outtxt
    errcounts = defaultdict(int)
    for l in alltxt.split('\n'):
        if 'err' in l:
            errcounts[l] += 1
    for k in errcounts.keys():
        out += f'{errcounts[k]}x: {k}\n'

    # Show full contents of .err and .out files
    out += f'\n{outstub}.err:\n\n{errtxt}'
    out += f'\n{outstub}.out:\n\n{outtxt}'
    return out
    
def exit_if_not_present(filename, email=False):
    if not os.path.exists(filename):
        print(f'Error: No {filename}. Exiting.')
        if email:
            dir = os.getcwd().split('SPARC/')[-1]
            body = f'Error: No {filename}\n\n'
            body += f'https://www1.psfc.mit.edu/~sballinger/solps/{path[-2]}/{path[-1]}\n\n'
            body += fail_info()
            send_email(f'Run stopped in {dir}', body)
        sys.exit(0)
        
def print_job_info(jobid):
    out = subprocess.getoutput(f'/bin/sacct -j {jobid} --format=JobID%9,JobName%45,Elapsed,State,NodeList%9')
    print('\n'.join(out.split('\n')[:3]))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('sims_to_run', type=int, help='Number of simulations to run.')
    parser.add_argument('-a', metavar='a', type=str, nargs='*', help='Annotation string describing any changes to the simulation for display in HTML index.')
    args = parser.parse_args()

    timequit = False
    userquit = False
    slownessquit = False
    notified_squeue_err = False
    refresh_seconds = 2*60
    
    if 'SOLPSTOP' not in os.environ.keys():
        print('Error: SOLPS not initialized. Exiting.')
        sys.exit(0)
    for filename in ['b2mn.dat','b2fstate','b2fstati']:
        exit_if_not_present(filename)
    
    signal.signal(signal.SIGINT, sigint_handler)
    signal.signal(signal.SIGQUIT, sigquit_handler)

    # Create run history directory if not present
    try:
        os.mkdir('history')
    except Exception:
        pass

    # Remove traces of prior failed run
    if os.path.exists('b2mn.exe.dir'):
        shutil.rmtree('b2mn.exe.dir')

    latest = get_latest_save()
    sims_to_run = int(args.sims_to_run)
    isim = 0

    if args.a != None:
        annotation = '\n'+str(latest+1) + ': ' + ' '.join(args.a)
        with open('history/annotations.txt','a') as f:
            f.write(annotation)

    for sim in range(latest+1, latest+sims_to_run+1):
        try:
            os.remove('current_run.txt')
        except:
            pass
        #subprocess.call(['quota', '-s', '-f', '/home'])
        path = os.getcwd().split('/')
        id = f'{path[-2]}_{path[-1]}_id{sim}'
        isim += 1
        print(f'Running {id} (run {isim}/{sims_to_run})')
        b2mn = parseDat('b2mn.dat')    
        total_iterations = int(b2mn["b2mndr_ntim"])
        # Not having b2fstate results in failure to launch
        subprocess.call(['cp', 'b2fstati', 'b2fstate'])
        process = subprocess.Popen(['psfcsubmit', '-m', '-np 32', '-j', id], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # Get SLURM job id
        allout = ''
        while 'Submitted batch job' not in allout:
            output = process.stdout.readline()
            if output:
                print(output.strip().decode('utf-8'))  # decode bytes and strip extra spaces
                allout += output.decode('utf-8')
            time.sleep(0.1)
        jobid = str(int(allout.split('Submitted batch job ')[-1]))
        with open('current_run.txt', 'w') as f:
            f.write(f'id{sim} {jobid}')
        # Could take a while for job to be accepted
        started_waiting_time = time.time()
        remaining_check_time = time.time()
        with Progress(TextColumn("[progress.description]{task.description}"), 
                      TimeElapsedColumn()) as progress:
            task = progress.add_task("Waiting for job to start", total=10, start=True)
            # Wait for job to start running. Exit if running or if case completed very quickly.
            # Keep looping if js is 'error' or ''
            js = job_status(jobid)
            while js != 'RUNNING' and js != 'COMPLETED':
                # Poll queue every 5 seconds in the first 5 min, then every 5 minutes to not overload
                if time.time()-started_waiting_time < 5*60:
                    time.sleep(5+random.randint(0,5))
                else:
                    time.sleep(5*60+random.randint(0,30))
                # Every 20 minutes, check est. start time and number of jobs ahead in queue
                if time.time()-remaining_check_time > 20*60:
                    try:
                        eststart = subprocess.getoutput(f'squeue -j {jobid} --start --format="%S" --noheader')
                        numahead = subprocess.getoutput(f'squeue_ahead {jobid}')
                        print(f'Estimated start time is {eststart}, {numahead} jobs ahead in queue')
                        remaining_check_time = time.time()
                    except Exception as e:
                        print('Error checking queue est. start and num. ahead:', e)
                js = job_status(jobid)
        time.sleep(10)
        print_job_info(jobid)
        # Wait for job to finish
        started = False
        timequit = False
        seconds_remaining = 0
        iterations_completed = 0
        start_time = 0
        print('Press ctrl-C to scancel job, ctrl-\ to quit cleanly and store')
        with Progress(TextColumn("[progress.description]{task.description}"),
                    BarColumn(bar_width=21),
                    MofNCompleteColumn(),
                    TaskProgressColumn(),
                    TimeElapsedColumn(),
                    TimeRemainingColumnCustom()) as progress:
            task = progress.add_task(f"id{sim} ({isim}/{sims_to_run})", total=total_iterations, start=False)
            refresh_seconds = 10
            
            # While job is running, track status. "not 'COMPLETED'" is resilient to errors checking status from sacct/slurm
            while job_status(jobid) != 'COMPLETED':
                try:
                    if os.path.isfile("b2mn.exe.dir/b2ftrace"):
                        iterations_completed = getIterations()
                        if iterations_completed > 0:
                            if not started:
                                progress.start_task(task)
                                started = True
                                start_time = time.time()
                                start_iterations = iterations_completed
                            # After start, check progress every 5 minutes plus random between 1 and 30 seconds to avoid traffic jams
                            if iterations_completed-start_iterations > 0:
                                refresh_seconds = 5*60+random.randint(0,30)
                                seconds_remaining = (time.time()-start_time)/(iterations_completed-start_iterations)*(total_iterations-iterations_completed)
                            # Exit cleanly and produce plots if near 6 hour job limit
                            if time.time()-start_time > 5.5*60*60 and not timequit:
                                subprocess.call(['touch','b2mn.exe.dir/.quit'])
                                timequit = True
                                refresh_seconds = 10
                                if iterations_completed < 100:
                                    slownessquit = True
                            last_iter_check_time = time.time()
                        progress.update(task, completed=iterations_completed)
                except Exception as e:
                    print('Track iterations exception:', e) 
                time.sleep(refresh_seconds)
            progress.update(task, completed=getIterations())
            # if not timequit and not userquit:
            #     progress.update(task, completed=total_iterations)
            
        # Print job info such as finish status, elapsed time
        print_job_info(jobid)
        # subprocess.call(['post', str(sim)])
        plotfile = f'history/{sim}/plots{sim}.pdf'
        # If making plots failed or b2fstate not present, something needs to be fixed
        exit_if_not_present(plotfile, email=True)
        exit_if_not_present('b2fstate', email=True)
        print(f'https://www1.psfc.mit.edu/~sballinger/solps/{path[-2]}/{path[-1]}')
        # If case blew up, stop and email me
        # try:
        #     send_email()
        # except Exception as e:
        #     print('')
        subprocess.call(['cp', 'b2fstate', 'b2fstati'])
        print()
        if userquit or slownessquit:
            break
        
if slownessquit:
    dir = os.getcwd().split('SPARC/')[-1]
    print('Run stopped due to slowness')
    send_email(f'Run stopped for slowness in {dir}', f'https://www1.psfc.mit.edu/~sballinger/solps/{path[-2]}/{path[-1]}')

if not userquit and not slownessquit:
    # Notify user that run is complete
    dir = os.getcwd().split('SPARC/')[-1]
    send_email(f'Run complete in {dir}', f'https://www1.psfc.mit.edu/~sballinger/solps/{path[-2]}/{path[-1]}')
