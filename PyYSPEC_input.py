#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  PyYSPEC_input.py
#   Purpose:   Automatically creates yspec input file for an event directory
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
except: sys.exit('usage: python PyYSPEC_input.py <path/to/yspec_infiles>')
indir = sys.argv[1]
create_source_inp(indir=indir)
