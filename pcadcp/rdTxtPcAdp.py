'''
Created on Nov 27, 2014

@author: bogdan
'''

import csv, os
import numpy
from datetime import datetime
from matplotlib.dates import date2num, num2date
import matplotlib.dates as dates
import re

class PcADP(object):
    '''
    classdocs
    '''



    def __init__(self, name, path, ctlname, filenames, num_segments, tinterv):
        '''
        Constructor
        '''
        self.profiles = {"start":"Time of first profile -------->",
              "end":"Time of last profile --------->",
              "dt":"Profile Interval ------------->" }

        self.path = path
        self.filenames = filenames
        self.num_segments = num_segments
        self.readCTL(ctlname)
        self.sitename = name
        self.timeinterv = tinterv
        self.goodbins=0

    def readCTL(self, ctlname):

        ctl = open(self.path + '/' + ctlname, "r")

        for line in ctl:
            if re.match(self.profiles["start"], line):
                startdate = line[line.index('>') + 2:]
            if re.match(self.profiles["end"], line):
                enddate = line[line.index('>') + 2:]
            if re.match(self.profiles["dt"], line):
                ends = line.index('s')
                dt = line[line.index('>') + 2:ends - 1]

        dtime = datetime.strptime(startdate.rstrip(), "%Y/%m/%d %H:%M:%S")
        start_num = dates.date2num(dtime)
        dtime = datetime.strptime(enddate.rstrip(), "%Y/%m/%d %H:%M:%S")
        end_num = dates.date2num(dtime)
        self.timespan = [start_num, end_num]
        self.sampletime = float(dt.rstrip()) / 3600. / 24.  # [s] must be converted to days




    def readVel(self):
        '''
         Sensor depth 0.5 m - distances below are measured from the sensor
         Column   Contents
         1        Profile #
         2        Velocity Cel 1  0.3-1.3m [units are in cm/s]
         3        Velocity Cel 2  1.3-2.3m [units are in cm/s]
         4        Velocity Cel 3  2.3-3.3m [units are in cm/s]
         5        Velocity Cel 4  3.3-4.3m [units are in cm/s]
         6        Velocity Cel 5  4.3-5.3m [units are in cm/s]

        '''

        if self.timeinterv != None:
            dtime = datetime.strptime(self.timeinterv[0], "%y/%m/%d %H:%M:%S")
            startt = dates.date2num(dtime)
            dtime = datetime.strptime(self.timeinterv[1], "%y/%m/%d %H:%M:%S")
            endt = dates.date2num(dtime)

        velocities = []
        timevec = []
        i = 0
        for filename in self.filenames:
            ifile = open(self.path + '/' + filename, 'rt')
            reader = csv.reader(ifile, skipinitialspace = True, delimiter = ' ', quotechar = '"')

            rownum = 0
            velocities.append([])

            time = self.timespan[0]

            for row in reader:
                if not row:
                    continue


                if self.timeinterv != None:
                    if time < startt or time > endt:
                            time += self.sampletime
                            continue

                velocities[i].append([])
                if i == 0:
                    timevec.append([])

                try:
                    colno = 0
                    for col in row:
                        if colno != 0:
                            # Convert:data is in c,/s and results are in m/s
                            velocities[i][rownum].append(float(col) / 100.0)
                            if i == 0:
                                timevec[rownum].append(time)
                        colno += 1

                    time += self.sampletime
                    rownum += 1
                except:
                    print("Error: read_data")
            ifile.close()
            i += 1
        # transpose
        tvelocities = {}
        j = 0
        for ar in range(0, len(velocities)):
            nar = numpy.array(velocities[ar], dtype = float)
            transp = numpy.transpose(nar)
            fname, ext = os.path.splitext(os.path.basename(self.filenames[j]))
            tvelocities[ext[1:]] = transp
            j += 1
        timearr = numpy.array(timevec, dtype = float)
        timearrt = numpy.transpose(timearr)
        self.goodbins= len(velocities[0][0])  # -1  # for this case thsi is not necessary
        return tvelocities, numpy.array(timearrt), self.goodbins 

    

    def readSTD(self):
        pass

    def readSpd(self):
        pass

    def readSNR(self):
        '''
        Signal to noise ratio
        '''
        pass

    def readTemp(self):
        [pn, Y, M, D, h, m, s, sampl, sndspd, head, pitch, roll, temp, pres, stdhd, stdpit, stdrol, stdtmp, stdpres, volt, time] = \
        self.readHDR()
        return [time, temp]

    def readHDR(self):
        '''
        Column    Contents                                                Units
        ------------------------------------------------------------------------------
        1         Profile number in file
        2         Profile time (start of averaging interval ) - Year
        3         Profile time (start of averaging interval ) - Month
        4         Profile time (start of averaging interval ) - Day
        5         Profile time (start of averaging interval ) - Hour
        6         Profile time (start of averaging interval ) - Minute
        7         Profile time (start of averaging interval ) - Second
        8         Number of samples averaged for this profile
        9         Sound speed to calculate velocity                        m/s
        10        Mean heading                                            degrees
        11        Mean pitch (rotation about the Y axis)                  degrees
        12        Mean roll (rotation about the X axis)                   degrees
        13        Mean temperature                                        deg C
        14        Mean pressure                                           dBar
        15        Standard deviation of heading                           degrees
        16        Standard deviation of pitch                             degrees
        17        Standard deviation of roll                              degrees
        18        Standard deviation of temperature                       deg C
        19        Standard deviation of pressure                          dBar
        20        Instrument power supply voltage level                   V
        ----------------------------------------------------------------------------
        21          ExtSensor1 (counts)
        22          ExtSensor2 (counts)
        23          ExtSensor1 Std Dev (counts)
        24          ExtSensor2 Std Dev (counts)

        '''
        if self.timeinterv != None:
            dtime = datetime.strptime(self.timeinterv[0], "%y/%m/%d %H:%M:%S")
            startt = dates.date2num(dtime)
            dtime = datetime.strptime(self.timeinterv[1], "%y/%m/%d %H:%M:%S")
            endt = dates.date2num(dtime)

        i = 0
        for filename in self.filenames:
            ifile = open(self.path + '/' + filename, 'rb')
            reader = csv.reader(ifile, skipinitialspace = True, delimiter = ' ', quotechar = '"')

            rownum = 0
            pn = []
            Y = [] 
            M= []
            D= []
            h= []
            m= []
            s= []
            sampl= []
            sndspd= []
            head= []
            pitch= []
            roll= []
            temp= []
            pres= []
            stdhd= []
            stdpit= []
            stdrol= []
            stdtmp= []
            stdpres= []
            volt= []
            time = []

            for row in reader:
                if not row:
                    continue

                try:
                    dtime = datetime(int(row[1]),int(row[2]), int(row[3]), int(row[4]), int(row[5]), int(row[6]))
                    dtnum = dates.date2num(dtime)
                    if self.timeinterv != None:
                        if dtnum < startt or dtnum > endt:
                            continue
                        
                    time.append(dtnum)
                    pn.append(int(row[0]))
                    Y.append(int(row[1])) 
                    M.append(int(row[2]))
                    D.append(int(row[3]))
                    h.append(int(row[4]))
                    m.append(int(row[5]))
                    s.append(int(row[6]))
                    sampl.append(int(row[7]))
                    sndspd.append(float(row[8]))
                    head.append(float(row[9]))
                    pitch.append(float(row[10]))
                    roll.append(float(row[11]))
                    temp.append(float(row[12]))
                    pres.append(float(row[13]))
                    stdhd.append(float(row[14]))
                    stdpit.append(float(row[15]))
                    stdrol.append(float(row[16]))
                    stdtmp.append(float(row[17]))
                    stdpres.append(float(row[18]))
                    volt.append(float(row[19]))
                except:
                    print("Error: read_data")
                
              
                
                rownum += 1
                
            ifile.close()
            i += 1
        return  [pn, Y, M, D, h, m, s, sampl, sndspd, head, pitch, roll, temp, pres, stdhd, stdpit, stdrol, stdtmp, stdpres, volt, time]

