#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


'/Users/vijays/Documents/environments/catmap_sigma/bin/python' < 'mkm_job.py' > 'aiida.out' 
