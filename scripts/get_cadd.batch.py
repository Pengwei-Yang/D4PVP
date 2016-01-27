#! /usr/bin/env python

import os
import sys
import subprocess
import datetime
import time
import signal

inFn = sys.argv[1]
outFn = sys.argv[2]

batchSize = 100

nLines = len(open(inFn, 'r').readlines())

nBatches = nLines / batchSize

leftOver = nLines - (batchSize * nBatches)
nBatches += 1 if leftOver > 0 else 0

print nBatches

headerIdx = 1
for batchIdx in range(0, nBatches):
    startIdx = batchSize * batchIdx
    endIdx = (batchSize * batchIdx) + batchSize

    if batchIdx == (nBatches - 1):
        endIdx = nLines
    if batchIdx > 0:
        headerIdx = 0
  
    batchOutFn = outFn + '.' + str(batchIdx)
   
    print batchIdx, startIdx, endIdx
    cmd = 'python get_cadd.py %s %s %s %s-%s' % (inFn, batchOutFn, str(headerIdx), str(startIdx), str(endIdx))
    print cmd

    os.system(cmd)
    print 'Command issued'
