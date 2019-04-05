#!/usr/bin/env python
import os
import sys
import popen2

prefices = dict([(x.split('_')[0],'') for x in popen2.popen3('ls')[0].read().strip().split('\n') if x != sys.argv[0]]).keys()

for prefix in prefices:
    cmd = 'cat ' + prefix + "* > " + prefix + ".fq"
    os.system(cmd)

