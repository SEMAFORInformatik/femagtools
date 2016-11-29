#!/usr/bin/env python
#
import femagtools.job
import tempfile
import os


def test_condor():
    workdir = tempfile.mkdtemp()
    job = femagtools.job.CondorJob(workdir)
    task = job.add_task()
    task.add_file('femag.fsl', ['exit_on_end=True'])
    try:
        job.prepareDescription()
    
        with open(os.path.join(workdir, 'femag.submit')) as f:
            x = [l.strip().split('=') for l in f]

        len(x) == 12
        x[-1][0].split() == '1'

        d = {l[0].strip(): l[1].strip() for l in x if len(l) > 1}
        d['Universe'] == 'vanilla'
    except Exception as e:
        pass # no femag installed on ci server.
    

def test_job():
    workdir = tempfile.mkdtemp()
    job = femagtools.job.Job(workdir)
    task = job.add_task()
    task.add_file('femag.fsl', ['exit_on_end=True'])
    
    with open(os.path.join(workdir, '0/femag.fsl')) as f:
        x = [l.strip().split('=') for l in f]
    d = {k: v for k, v in x}
    d['exit_on_end'] == 'True'
