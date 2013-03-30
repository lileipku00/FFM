from MT_handler import event_finder
from obspy.core import UTCDateTime
t = UTCDateTime('2009-09-30-10-16-09')
ev = event_finder(t, -0.72, 99.87, 7.5)
ev.search_neic('../EventCatalogues/mt.neic.bin')
