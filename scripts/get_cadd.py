#! /usr/bin/env python

import os
import sys
import subprocess
import datetime
import time
import signal

clinvarFn = sys.argv[1]
clinvarOutFn = sys.argv[2]
addHeaderParam = sys.argv[3]
lineRange = sys.argv[4]

caddHttp = 'http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs_inclAnno.tsv.gz'

def timeout_command(command, timeout):
    """call shell-command and either return its output or kill it
    if it doesn't normally exit within timeout seconds and return None"""

    start = datetime.datetime.now()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    while process.poll() is None:
        time.sleep(0.1)
        now = datetime.datetime.now()
        if (now - start).seconds > timeout:
            os.kill(process.pid, signal.SIGKILL)
            os.waitpid(-1, os.WNOHANG)
            return None
    return process.stdout.read()

clinvarFlines = open(clinvarFn, 'r').readlines()
clinvarOutFile = open(clinvarOutFn, 'w')

addHeader = False
if addHeaderParam == '1':
    addHeader = True

start = 0
end = len(clinvarFlines)

if lineRange.find('-') > -1:
    lRange = lineRange.split('-')
    start = int(lRange[0])
    end = int(lRange[1])

idx = 0
for line in clinvarFlines[start:end]:
    if line.find('#') > -1:
        continue
    line = line.strip()
    linesplit = line.split('\t')

    chr = linesplit[0]
    pos = linesplit[1]
    
    posStr = chr + ':' + str(pos) + '-' + str(pos)
    if addHeader:
        cmd = 'tabix -h %s %s' %(caddHttp, posStr)
        addHeader = False
    else:
        cmd = 'tabix %s %s' % (caddHttp, posStr)

    success = False
    nFails = 0
    while not success:
        output = timeout_command(cmd, 2)
        if output is not None:
            success = True
        else:
            nFails += 1
            print nFails, 'fails'
#    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
#    outStr, outErr = p.communicate()
#    print outStr, outErr
    clinvarOutFile.write(output)
    print idx, posStr
    idx += 1
