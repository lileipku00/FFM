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

# Required Python and Obspy modules will be imported in this part.

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import sys

from util_PyYSPEC import *


try: print 'input:\n%s' %(sys.argv[1])
except: sys.exit('usage: python PyYSPEC <path/to/yspec_infiles>')
indir = sys.argv[1]
create_source_inp(indir=indir)

try: print 'input:\n%s %s %s %s' %(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
except: sys.exit('ERROR in the number of inputs!')
run_yspec(indir_submit=sys.argv[2], yspec_inp=sys.argv[3], 
            num_proc=sys.argv[4], indir_output=sys.argv[5])

