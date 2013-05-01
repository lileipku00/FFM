#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  PyYSPEC.py
#   Purpose:   Automatically runs YSPEC for an event directory
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#
#   Copyright (C) 2013 Kasra Hosseini
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import filecmp
import sys, os
import shutil

from util_PyYSPEC import *

'''
sys.argv[1]: path to yspec infiles/yspec directory
sys.argv[2]: path to bin directory of yspec
sys.argv[3]: path to yspec input file (yspec.in) in infiles/yspec
sys.argv[4]: number of processes
sys.argv[5]: path to YSPEC gallery (YSPEC_SYN_GALLERY)
'''

try: print 'input:\n%s' %(sys.argv[1])
except: sys.exit('usage: python PyYSPEC <path/to/yspec_infiles>')
indir = sys.argv[1]
create_source_inp(indir=indir)

try: print 'input:\n%s %s %s %s' %(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
except: sys.exit('ERROR in the number of inputs!')

if os.path.isfile(os.path.join(sys.argv[5], 'yspec.in')):
    if not filecmp.cmp(sys.argv[3], os.path.join(sys.argv[5], 'yspec.in')):
        print '\nCurrent event was simulated before with different setting.'
        print 'Removing the previous one and start the new simulation.'
        print '\nRemoved directory:'
        print sys.argv[5]
        shutil.rmtree(sys.argv[5])
    else:
        print '\nThe simulation was done before and it is found in the archive:'
        print sys.argv[5]
else:
    print '\nStart a new simulation for:'
    print sys.argv[5]
run_yspec(indir_submit=sys.argv[2], yspec_inp=sys.argv[3], 
            num_proc=sys.argv[4], indir_output=sys.argv[5])

