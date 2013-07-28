band = '01'
minxcor = 0.85
lonang = 2
latang = 1
adrsdir = '/import/neptun-helles/hosseini/FFM/'

#==================================================
# Function GETEVLOC to get the location of the event
#===================================================
def getevloc(event):
        id = open(event + '/outfiles/ffproc.source')
        lines = id.readlines()
        temp = lines[3]
        temp = temp.split()
        evlat = float(temp[0])
        evlon = float(temp[1])
        
        return evlat, evlon
#====================================================
# Function Reads in the station DATA, for those > xcor
#====================================================
#
#       1- compares the xcor of each station with xmincor
#       2- and appends [dt,statlat, statlon] to eventdata
#               retunrs:        eventdata[dt, statlat, statlon]
def getstat(event = '/import/neptun-radler/hosseini-downloads/KASRA/FFM/0274.2005.100.a', band = '01', minxcor = 0.9):
        id = open(event + '/outfiles/ffproc.ampstt.band'+band, 'r')
        line = id.readlines()
        eventdata = []
        for i in range(2, len(line)):
                temp = line[i]
                temp = temp.split()
                xcor = float(temp[2])                   # step 1
                if xcor < minxcor:
                        continue
                dt = float(temp[5])
                statlat = float(temp[6])                # step2
                statlon = float(temp[7])
                eventdata.append([dt, statlat, statlon])
        return eventdata
#===================================================
# Function Midpoint
# calculates the midpoint based on spherical geometry
def midpoint(stlat, stlon, evlat, evlon):
        from numpy import deg2rad, rad2deg, arctan2, sin, cos, sqrt
        evlon = deg2rad(evlon)
        evlat = deg2rad(evlat)
        stlon = deg2rad(stlon)
        stlat = deg2rad(stlat)
        bx = cos(stlat)*cos(stlon-evlon)
        by = cos(stlat)*sin(stlon-evlon)
        midlat = arctan2(sin(evlat)+sin(stlat),sqrt((cos(evlat)+bx)*(cos(evlat)+bx)+by*by) )
        midlat = rad2deg(midlat)
        if midlat >=90.00:
                midlat = midlat - 180.00
        elif midlat <= -90.00:
                midlat = midlat + 180.00
        midlon = evlon + arctan2(by , cos(evlat)+bx)
        midlon = rad2deg(midlon)
        if midlon >=180.00:
                midlon = midlon - 360.00
        elif midlon <= -180.00:
                midlon = midlon + 360.00
        
        return midlat, midlon
#====================================================
# function                      STATDATA:(event, band)
#                                       return :        eventdata, element
#       
#       1- event's lat & lon = getevloc(event)
#       2- station's data : eventdata = getstat(event, band, minxcor)
#       3- calculate the midpoint and replace eventdata[][1] & eventdata[][2]
#       4- subtract the calculated midpoint from all dt for this event
#       5- gets the correspoding address of the data:   element = elemdist(eventdata, lonang, latang)
def statdata(event, band = '01'):
        import numpy as np
        import math

        try:
                eventdata = getstat(event, band, minxcor)                       # step 2
                evlat, evlon = getevloc(event)                                 # step 1
                dt = []
                for ii in range(len(eventdata)):                                        # step 3
                        eventdata[ii][1],eventdata[ii][2]=midpoint(eventdata[ii][1],eventdata[ii][2], evlat, evlon)
                        dt.append(eventdata[ii][0])
                if np.mod(len(dt),2) == 1:                                              #step 4
                        dtmed = np.median(dt)
                else:
                        dt.append(dt[-1])
                        dtmed = np.median(dt)
                for jj in range(len(eventdata)):
                        eventdata[jj][0] = eventdata[jj][0] - dtmed
                element = elemdist(eventdata, lonang, latang)                           #step 5

                return eventdata, element
        except Exception, e:
                print e
#===================================================
# function elemdist, finds the element for each midpoint of stations
#       1 - define an int to shift all cells to positive addresses
#       2 - devide by the input angles, take the integer,
#                       if addresses were positive shift FWD with latconst
#                       if negative subtract 1, then shift FWD          ---> to prevent overlap
#                                                                       of first positive and negative cells
def elemdist(eventdata, lonang = 0.5, latang = 0.5 ):
        latconst = int(90/latang)
        lonconst = int(180/lonang)
        element = []
        for kk in range(len(eventdata)):
                if eventdata[kk][1] >= 0:
                        latel = int(eventdata[kk][1]/latang)+latconst
                        if eventdata[kk][2] >= 0:
                                lonel = int(eventdata[kk][2]/lonang)+lonconst
                        else:
                                lonel = int(eventdata[kk][2]/lonang)-1 + lonconst
                else:
                        latel = int(eventdata[kk][1]/latang)-1+latconst
                        if eventdata[kk][2] >= 0:
                                lonel = int(eventdata[kk][2]/lonang)+lonconst

                        else:
                                lonel = int(eventdata[kk][2]/lonang)-1 + lonconst
                
                temp = '0'*(3-len(str(latel))) + str(latel)+ '0'*(3-len(str(lonel)))+str(lonel)
                element.append(temp)
        return element
#===================================================
# builds a matrix with the dimensions calculated
#def matrixbuild(lonang = 60, latang = 30, eventdata, element)
#        maxlong = 180/lonang
#        minlong = -180/lonang
#        maxlat = 90/latang
#        minlat = -90/latang
#main program
import numpy as np
import glob as gl
event = gl.glob(adrsdir+'0*')
eventdata = []
elements = []
all_event_element = []
for jj in range(len(event)):
        try:
                eventdatatemp, elementstemp = statdata(event[jj], band)
                for bb in range(len(eventdatatemp)):
                        eventdata.append(eventdatatemp[bb])
                        elements.append(elementstemp[bb])
                        all_event_element.append([eventdatatemp[bb],elementstemp[bb]])
        except Exception, e:
                print e
                continue

enum = 0
all_event_element.sort(key=lambda x: x[1])
for i in range(len(all_event_element)-1, -1, -1):
        if all_event_element[i][1] == all_event_element[i-1][1]:
                all_event_element[i-1][0][0] = all_event_element[i][0][0] + all_event_element[i-1][0][0]
                del all_event_element[i]
                enum += 1
        else:
                if enum > 0:
                        all_event_element[i][0][0] = all_event_element[i][0][0]/enum
                        enum = 0
final_list = []
latconst = int(90/latang)
lonconst = int(180 / lonang)
for i in range(len(all_event_element)):
        dt = all_event_element[i][0][0]
        elem = all_event_element[i][1]
        lat = (int(elem[0:3])-latconst)*latang + (latang/2)
        lon = (int(elem[3:])-lonconst)*lonang + (lonang/2)

        final_list.append([lat,lon, dt])

#fout = open('result', 'w')
#fout.write('latitude  lontitude  dT\n')
#for jj in range(len(final_list)):
#        #temp = ' '*(8-len(str(final_list[jj][0]))) + str(final_list[jj][0]) + \
#        #        ' '*(8-len(str(final_list[jj][1])))+str(final_list[jj][1]) + \
#        #        ' '*(8-len(str(final_list[jj][2])))+str(final_list[jj][2])+'\n'
#        temp = '%s, %s, %s\n' %(final_list[jj][0], final_list[jj][1], final_list[jj][2])
#        fout.write(temp)
#fout.close()

print"""
def midpoint(stlat, stlon, evlat, evlon):
        return: midlat, midlon
def getevloc(event):
        return: evlat, evlon
def getstat(event = '/import/neptun-radler/hosseini-downloads/KASRA/FFM/0274.2005.100.a', band = '01', minxcor = 0.9):
        return: eventdata.append([dt, statlat, statlon])
def elemdist(eventdata, lonang = 0.5, latang = 0.5 ):
        return element
all_event_element[eventdata,element]
final_list[lat,lon,dt]
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

m = Basemap(projection='moll', lon_0=-90, lat_0=0, resolution='c')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary()

map_x = []
map_y = []
map_dt = []
for i in range(len(final_list)):
    if -5<=final_list[i][2]<=5:
        x, y = m(final_list[i][1], final_list[i][0])
        map_x.append(x)
        map_y.append(y)
        map_dt.append(final_list[i][2])
m.scatter(map_x, map_y, c=map_dt,vmin=min(map_dt), vmax=max(map_dt),alpha=0.6, edgecolor='none', zorder=10)
m.colorbar()
plt.show()

