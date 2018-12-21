#!/usr/bin/env python3
"""
Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --cluster slurm_scheduler.py --cluster-status slurm_status.py
"""

import os
import re
import sys
import warnings
import subprocess

i = 1
while not re.match("\d+", sys.argv[i]):
  i=i+1
jobid = sys.argv[i]


out= subprocess.run(['scontrol','show','jobid',jobid],stdout=subprocess.PIPE).stdout.decode('utf-8')

def parse_key_value(stream):
    params={}
    for key_value_pair in stream.split():
        name, var = key_value_pair.partition("=")[::2]
        params[name.strip()] = var
    return params

state = "FAILED"

if 'JobState' in parse_key_value(out).keys():
  state = parse_key_value(out)['JobState']

map_state={"PENDING":'running',
           "RUNNING":'running',
           "SUSPENDED":'running',
           "CANCELLED":'failed',
           "COMPLETING":'running',
           "COMPLETED":'success',
           "CONFIGURING":'running',
           "FAILED":'failed',
           "TIMEOUT":'failed',
           "PREEMPTED":'failed',
           "NODE_FAIL":'failed',
           "REVOKED":'failed',
           "SPECIAL_EXIT":'failed',
           "":'success'}

print(map_state[state])