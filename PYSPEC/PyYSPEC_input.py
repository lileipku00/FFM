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
import sys

from util_PyYSPEC import *

"""
sys.argv[1]: path to yspec infiles/yspec directory
"""

try:
    print 'input:\n%s' % sys.argv[1]
except Exception, e:
    sys.exit('usage: python PyYSPEC_input.py <path/to/yspec_infiles>\nERROR: %s' % e)

indir = sys.argv[1]
create_source_inp(indir=indir)
