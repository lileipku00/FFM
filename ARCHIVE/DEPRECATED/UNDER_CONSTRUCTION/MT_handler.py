#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  MT_handler.py
#   Purpose:   Read NEIC and Harvard MT catalog
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#
#   Copyright (C) 2013 Kasra Hosseini
#-------------------------------------------------------------------

# Added this line for python 2.5 compatibility
from __future__ import with_statement
from obspy.core import UTCDateTime

class event_finder:
    def __init__(self, ev_datetime, ev_lat, ev_lon, ev_mag):
        self.evref_utc=ev_datetime
        self.evref_lat=ev_lat
        self.evref_lon=ev_lon
        self.evref_mag=ev_mag
    def search_neic(self, neicpath):
        neic_open=open(neicpath)
        self.neic_read=neic_open.readlines()
        for i in xrange(len(self.neic_read)):
            self.neic_read[i] = self.neic_read[i].split()
        found_event_0=[]
        for i in xrange(len(self.neic_read)):
            if self.evref_utc.year==int(self.neic_read[i][0]):
                if self.evref_utc.month == int(self.neic_read[i][1]):
                    neic_utc=UTCDateTime(self.neic_read[i][0]+'-'+self.neic_read[i][1]+'-'+\
                                         self.neic_read[i][2]+'T'+self.neic_read[i][3])
                    if abs(self.evref_utc-neic_utc) < 60.0:
                        found_event_0.append([i, self.neic_read[i]])
        if len(found_event_0)==1:
            print 'Event is found with the following information:'
            print found_event_0[0][1]
            found_event_1 = found_event_0
        elif len(found_event_0)==0:
            print 'search_neic could not find any event! try search_harvard!'
            found_event_1=[]
        elif len(found_event_0)>1:
            found_event_1=[]
            print '%s events found based on the datetime.' %len(found_event_0)
            for i in xrange(len(found_event_0)):
                if abs(self.evref_mag-float(self.neic_read[found_event_0[i][0]][8])) < 0.2:
                    found_event_1.append(found_event_0[i][0],self.neic_read[found_event_0[i][0]])
        self.neic_events=found_event_1
