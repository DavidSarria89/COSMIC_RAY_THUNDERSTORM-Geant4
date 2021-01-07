import platform
import subprocess as sp
from subprocess import call
import time
from os import system
from os import walk
from os import listdir
from mpi4py import MPI
from functools import partial
from multiprocessing.dummy import Pool
import sys

computer_name = platform.node()

IS_CLUSTER = False  # use of MPI

################################################################################
# functions definitions for MPI usage (if Cluster)


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

################################################################################
# Defining commands to run


nb_run = 100
POTENTIAL_LIST = [    0,    25,    50,  100,  150,  200]
STATS_LIST =     [10000,  2500,  2500,   20,   20,   20]

# defining the commands to be run in parallel
commands = []
excecutable = './mos_test'

for _ in range(nb_run):
    for ii, POT in enumerate(POTENTIAL_LIST):
        commands.append(excecutable + ' ' + str(POT) + ' ' + str(STATS_LIST[ii]))

################################################################################
# LOCAL RUN (uses python multiprocessing library)

if not IS_CLUSTER:
    nb_thread = 4  # number of threads (cpu) to use

    # Making an array where each element is the list of command for a given thread

    command_number = len(commands)

    print('Number of commands required ' + str(command_number))

    pool = Pool(nb_thread)  # to be always set to 1 for this MPI case
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))

# MPI run (uses mpi4py)
else:
    # MPI initializations and preliminaries
    comm = MPI.COMM_WORLD   # get MPI communicator object
    size = comm.Get_size()       # total number of processes
    rank = comm.Get_rank()       # rank of this process
    status = MPI.Status()   # get MPI status object
    # automatically uses all CPU threads available

    # Define MPI message tags
    tags = enum('READY', 'DONE', 'EXIT', 'START')

    if rank == 0:
        # Master process executes code below
        tasks = commands
        task_index = 0
        num_workers = size - 1
        closed_workers = 0
        print("Master starting with %d workers" % num_workers)
        while closed_workers < num_workers:
            data = comm.recv(source=MPI.ANY_SOURCE,
                             tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            if tag == tags.READY:
                # Worker is ready, so send it a task
                if task_index < len(tasks):
                    comm.send(tasks[task_index], dest=source, tag=tags.START)
                    print("Sending task %d to worker %d" %
                          (task_index, source))
                    task_index += 1
                else:
                    comm.send(None, dest=source, tag=tags.EXIT)
            elif tag == tags.DONE:
                results = data
                print("Got data from worker %d" % source)
            elif tag == tags.EXIT:
                print("Worker %d exited." % source)
                closed_workers += 1

        print("Master finishing")
    else:
        # Worker processes execute code below
        name = MPI.Get_processor_name()
        print("I am a worker with rank %d on %s." % (rank, name))
        while True:
            comm.send(None, dest=0, tag=tags.READY)
            task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()

            if tag == tags.START:
                # Do the work here
                task2 = [task]
                pool = Pool(1)  # to be always set to 1 for this MPI case
                for i, returncode in enumerate(pool.imap(partial(call, shell=True), task2)):
                    if returncode != 0:
                        print("%d command failed: %d" % (i, returncode))
                comm.send(returncode, dest=0, tag=tags.DONE)
            elif tag == tags.EXIT:
                break

    comm.send(None, dest=0, tag=tags.EXIT)
