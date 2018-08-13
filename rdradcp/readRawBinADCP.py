'''
Code modified and speed improvements Bogdan Hlevca 14 November 2011 University of Toronto

def  readRawBinADCP(name, *varargin):
        returns [adcp,cfg,ens,hdr]
        Read (raw binary) RDI ADCP files,
        ADCP=RDRADCP(NAME) reads the raw binary RDI BB/Workhorse ADCP file NAME and
        puts all the relevant configuration and measured data into a data structure
        ADCP (which is self-explanatory). This program is designed for handling data
        recorded by moored instruments (primarily Workhorse-type but can also read
        Broadband) and then downloaded post-deployment. For vessel-mount data I
        usually make p-files (which integrate nav info and do coordinate transformations)
        and then use RDPADCP.

          This current version does have some handling of VMDAS, WINRIVER, and WINRIVER2 output
          files, but it is still 'beta'. There are (inadequately documented) timestamps
         of various kinds from VMDAS, for example, and caveat emptor on WINRIVER2 NMEA data.

          [ADCP,CFG]=RDRADCP(...) returns configuration data in a
          separate data structure.

          Various options can be specified on input:
          [..]=RDRADCP(NAME,NUMAV) averages NUMAV ensembles together in the result.
          [..]=RDRADCP(NAME,NUMAV,NENS) reads only NENS ensembles (-1 for all).
          [..]=RDRADCP(NAME,NUMAV,[NFIRST NEND]) reads only the specified range
           of ensembles. This is useful if you want to get rid of bad data before/after
           the deployment period.

          Notes- sometimes the ends of files are filled with garbage. In this case you may
                 have to rerun things explicitly specifying how many records to read (or the
                 last record to read). I don't handle bad data very well. Also - in Aug/2007
                 I discovered that WINRIVER-2 files can have a varying number of bytes per
                 ensemble. Thus the estimated number of ensembles in a file (based on the
                 length of the first ensemble and file size) can be too high or too low.

               - I don't read in absolutely every parameter stored in the binaries;
                 just the ones that are 'most' useful. Look through the code if
                 you want to get other things.

               - chaining of files does not occur (i.e. read .000, .001, etc.). Sometimes
                 a ping is split between the end of one file and the beginning of another.
                 The only way to get this data is to concatentate the files, using
                   cat file1.000 file1.001 > file1   (unix)
                   copy file1.000/B+file2.001/B file3.000/B     (DOS/Windows)

                 (as of Dec 2005 we can probably read a .001 file)

               - velocity fields are always called east/north/vertical/error for all
                 coordinate systems even though they should be treated as
                 1/2/3/4 in beam coordinates etc.

          String parameter/option pairs can be added after these initial parameters:
          'info':   'yes' or 'no' (dafault). yes will print additional time information and it is just a little slower

          'baseyear'    : Base century for BB/v8WH firmware (default to 2000).

          'despike'    : [ 'no' | 'yes' | 3-element vector ]
                         Controls ensemble averaging. With 'no' a simple mean is used
                         (default). With 'yes' a mean is applied to all values that fall
                         within a window around the median (giving some outlier rejection).
                         This is useful for noisy data. Window sizes are [.3 .3 .3] m/s
                         for [ horiz_vel vert_vel error_vel ] values. If you want to
                         change these values, set 'despike' to the 3-element vector.

         R. Pawlowicz (rich@eos.ubc.ca) - 17/09/99

         R. Pawlowicz - 17/Oct/99
                  5/july/00 - handled byte offsets (and mysterious 'extra" bytes) slightly better, Y2K
                  5/Oct/00 - bug fix - size of ens stayed 2 when NUMAV==1 due to initialization,
                             hopefully this is now fixed.
                  10/Mar/02 - #bytes per record changes mysteriously,
                              tried a more robust workaround. Guess that we have an extra
                              2 bytes if the record length is even?
                  28/Mar/02 - added more firmware-dependent changes to format; hopefully this
                              works for everything now (put previous changes on firmer footing?)
                  30/Mar/02 - made cfg output more intuitive by decoding things.
                            - An early version of WAVESMON and PARSE which split out this
                              data from a wave recorder inserted an extra two bytes per record.
                              I have removed the code to handle this but if you need it see line 509
                 29/Nov/02  - A change in the bottom-track block for version 4.05 (very old!).
                 29/Jan/03  - Status block in v4.25 150khzBB two bytes short?
                 14/Oct/03  - Added code to at least 'ignore' WinRiver GPS blocks.
                 11/Nov/03  - VMDAS navigation block, added hooks to output
                              navigation data.
                 26/Mar/04  - better decoding of nav blocks
                            - better handling of weird bytes at beginning and end of file
                              (code fixes due to Matt Drennan).
                 25/Aug/04  - fixes to "junk bytes" handling.
                 27/Jan/05  - even more fixed to junk byte handling (move 1 byte at a time rather than
                              two for odd lengths.
          29/Sep/2005 - median windowing done slightly incorrectly in a way which biases
                        results in a negative way in data is *very* noisy. Now fixed.

           28/Dc/2005  - redid code for recovering from ensembles that mysteriously change length, added
                         'checkheader' to make a complete check of ensembles.
             Feb/2006  - handling of firmware version 9 (navigator)
           23/Aug/2006 - more firmware updates (16.27)
           23/Aug2006  - ouput some bt QC stiff
           29/Oct/2006 - winriver bottom track block had errors in it - now fixed.
           30/Oct/2006 - pitch_std, roll_std now uint8 and not int8 (thanks Felipe pimenta)
           13/Aug/2007 - added Rio Grande (firmware v 10),
                         better handling of those cursed winriver ASCII NMEA blocks whose
                         lengths change unpredictably.
                         skipping the inadequately documented 2022 WINRIVER-2 NMEA block
           13/Mar/2010 - firmware version 50 for WH.
           16/Nov/2011 - Translation to Numerical Python. numpy and scipy packages are required. mathplotlib the graphic package
                         There are also speed improvements (Bogdan Hlevca)
           02/Jul/2012-  Minor bug fixes and error handling (Bogdan HLevca)

====================================================================================================
NOTE: east velocity (u) is the east->west velocity;  North velocity (v) is the North->South velocity
====================================================================================================

'''

from . import scanf  # local file
import os, io, math
import numpy
import scipy

from matplotlib.dates import date2num
from datetime import datetime
import time


# globals - abd practice
hdr = None
ens = None
FIXOFFSET = None
SOURCE = None
cfg = None

class ADCPData:
    '''Place holder for the ADCP data
    '''

    def __init__(self, n, cfg):
        self.type = cfg.sourceprog

        self.goodbins = 0
        self.num_rec = 0
        self.name = 'adcp'
        self.config = cfg
        self.mtime = numpy.zeros((1, n))
        self.number = numpy.zeros((1, n))
        self.pitch = numpy.zeros((1, n))
        self.roll = numpy.zeros((1, n))
        self.heading = numpy.zeros((1, n))
        self.pitch_std = numpy.zeros((1, n))
        self.roll_std = numpy.zeros((1, n))
        self.heading_std = numpy.zeros((1, n))
        self.depth = numpy.zeros((1, n))
        self.temperature = numpy.zeros((1, n))
        self.salinity = numpy.zeros((1, n))
        self.pressure = numpy.zeros((1, n))
        self.pressure_std = numpy.zeros((1, n))
        self.east_vel = numpy.zeros((cfg.n_cells, n))
        self.north_vel = numpy.zeros((cfg.n_cells, n))
        self.vert_vel = numpy.zeros((cfg.n_cells, n))
        self.error_vel = numpy.zeros((cfg.n_cells, n))
        self.corr = numpy.zeros((cfg.n_cells, 4, n))
        self.status = numpy.zeros((cfg.n_cells, 4, n))
        self.intens = numpy.zeros((cfg.n_cells, 4, n))
        self.perc_good = numpy.zeros((cfg.n_cells, 4, n))
        self.bt_range = numpy.zeros((4, n))
        self.bt_vel = numpy.zeros((4, n))
        self.bt_corr = numpy.zeros((4, n))
        self.bt_ampl = numpy.zeros((4, n))
        self.bt_perc_good = numpy.zeros((4, n))

        if self.type == 'WINRIVER':
            self.nav_mtime = numpy.zeros((1, n))
            self.nav_longitude = numpy.zeros((1, n))
            self.nav_latitude = numpy.zeros((1, n))
        elif self.type == 'VMDAS':
            self.nav_smtime = numpy.zeros((1, n))
            self.nav_emtime = numpy.zeros((1, n))
            self.nav_slongitude = numpy.zeros((1, n))
            self.nav_elongitude = numpy.zeros((1, n))
            self.nav_slatitude = numpy.zeros((1, n))
            self.nav_elatitude = numpy.zeros((1, n))
            self.nav_mtime = numpy.zeros((1, n))


class ADCPCfg:
    ''' Place holder for the ADCP data configuration'''

    def __init__(self):
        self.name = ''
        self.sourceprog = 'instrument'  # default - depending on what data blocks are
        '''
         8,9,16 - WH navigator
         10 -rio grande
         15, 17 - NB
         19 - REMUS, or customer specific
         11- H-ADCP
         31 - Streampro
         34 - NEMO
         50 - WH, no bottom track (built on 16.31)
         51 - WH, w/ bottom track
         52 - WH, mariner
        '''
        self.prog_ver = ''
        self.config = None
        self.beam_angle = None
        self.numbeams = None
        self.beam_freq = None
        self.beam_pattern = None
        self.orientation = None
        self.simflag = None
        self.n_beams = None
        self.n_cells = None
        self.pings_per_ensemble = None
        self.cell_size = None
        self.blank = None
        self.prof_mode = None
        self.corr_threshold = None
        self.n_codereps = None
        self.min_pgood = None
        self.evel_threshold = None
        self.time_between_ping_groups = None
        self.coord = None
        self.coord_sys = None
        self.use_pitchroll = None
        self.use_3beam = None
        self.bin_mapping = None
        self.xducer_misalign = None
        self.magnetic_var = None
        self.sensors_src = None
        self.sensors_avail = None
        self.bin1_dist = None
        self.xmit_pulse = None
        self.water_ref_cells = None
        self.fls_target_threshold = None
        self.xmit_lag = None
        self.serialnum = None
        self.sysbandwidth = None
        self.syspower = None
        self.navigator_basefreqindex = None
        self.remus_serialnum = None
        self.h_adcp_beam_angle = None
        self.ranges = None
    # end

class Ensamble:
    def __init__(self, n, cfg):
        '''
            zeros with no classname/dtype specified is
            filled with numpy.float64 for matlab is not specified.
        '''
        self.number = numpy.zeros((1, n))
        self.rtc = numpy.zeros((7, n))
        self.BIT = numpy.zeros((1, n))
        self.ssp = numpy.zeros((1, n))
        self.depth = numpy.zeros((1, n))
        self.pitch = numpy.zeros((1, n))
        self.roll = numpy.zeros((1, n))
        self.heading = numpy.zeros((1, n))
        self.temperature = numpy.zeros((1, n))
        self.salinity = numpy.zeros((1, n))
        self.mpt = numpy.zeros((1, n))
        self.heading_std = numpy.zeros((1, n))
        self.pitch_std = numpy.zeros((1, n))
        self.roll_std = numpy.zeros((1, n))
        self.adc = numpy.zeros((8, n))
        self.error_status_wd = numpy.zeros((1, n))
        self.pressure = numpy.zeros((1, n))
        self.pressure_std = numpy.zeros((1, n))
        self.east_vel = numpy.zeros((cfg.n_cells, n))
        self.north_vel = numpy.zeros((cfg.n_cells, n))
        self.vert_vel = numpy.zeros((cfg.n_cells, n))
        self.error_vel = numpy.zeros((cfg.n_cells, n))
        # self.intens = numpy.zeros((n , cfg.n_cells, 4))
        # self.percent = numpy.zeros((n , cfg.n_cells, 4))
        # self.corr = numpy.zeros((n , cfg.n_cells, 4))
        # self.status = numpy.zeros((n , cfg.n_cells, 4))
        self.intens = numpy.zeros((cfg.n_cells, 4, n))
        self.percent = numpy.zeros((cfg.n_cells, 4, n))
        self.corr = numpy.zeros((cfg.n_cells, 4, n))
        self.status = numpy.zeros((cfg.n_cells, 4, n))
        self.bt_range = numpy.zeros((4, n))                 #BT = bottom track
        self.bt_vel = numpy.zeros((4, n))
        self.bt_corr = numpy.zeros((4, n))
        self.bt_ampl = numpy.zeros((4, n))
        self.bt_perc_good = numpy.zeros((4, n))
        self.smtime = numpy.zeros((1, n))
        self.emtime = numpy.zeros((1, n))
        self.slatitude = numpy.zeros((1, n))
        self.slongitude = numpy.zeros((1, n))
        self.elatitude = numpy.zeros((1, n))
        self.elongitude = numpy.zeros((1, n))
        self.nmtime = numpy.zeros((1, n))
        self.flags = numpy.zeros((1, n))

class Hdr:
    def __init(self):
        self.nbyte = numpy.nan
        self.dat_offsets = numpy.nan


# seek function instructions
bof = 0; cof = 1; eof = 2

def dec2bin(x):
    return x and (dec2bin(x / 2) + str(x % 2)) or ''
    # dec2bin = lambda x: (dec2bin(x/2) + str(x%2)) if x else ''

def dec2base(x, b):
    return x and ((dec2base(x / b, b, digits) + str(x % b)) or "")


######################################################################
# Decimal to any base conversion.
# Convert 'num' to a list of 'l' numbers representing 'num'
# to base 'base' (most significant symbol first).
######################################################################
def dec2base(num, base, l):
    s = ""
    n = num
    for i in range(l):
        s += str(n % base)
        n = int(n / base)
    if n != 0:
        print('Number ', num, ' requires more than ', l, 'digits.')
    # reverse the string because is bas build from low-> high bits
    return s[::-1]

######################################################################
# Conversion from any base to decimal.
# Convert a list 's' of symbols to a decimal number
# (most significant symbol first)
######################################################################
def base2dec(s, base):
    num = 0
    for i in range(len(s)):
        num = num * base + s[i]
    return num

# @countcalls
# def dec2base(x, b, digits):
#    if (dec2base.count() < digits) or y == "":
#        y = dec2base(x / b, b, digits) + str(x % b)
#    return y

def dec2hex(n):
    '''return the hexadecimal string representation of integer n'''
    return "%X" % n

def hex2dec(s):
    '''return the integer value of a hexadecimal string s'''
    return int(s, 16)

def strfind(longstr, substr):
    '''find all occurrences of a substring
    Returns an array of indexes wher the substr was found or empty array in no position was found
    '''
    return [m.start() for m in re.finditer(substr, longstr)]
#----------------------------------------
def checkheader(fd):
    '''
        Checks the header of the ADCP file
        returns boolean valid
    '''

    # disp('checking');
    valid = 0;

    numbytes = numpy.fromfile(fd, numpy.int16, 1)[0]  # Following the header bytes is numbytes

    if numbytes > 0:  # and we move forward numbytes>0
        try:
            fd.seek(numbytes - 2, cof)  # (numbytes-2,cof)
            cfgid = numpy.fromfile(fd, numpy.uint8, 2);
            if len(cfgid) == 2:  # will Skip the last ensemble (sloppy code)
                fd.seek(-numbytes - 2, cof)
                # # fprintf([dec2hex(cfgid(1)) ' ' dec2hex(cfgid(2)) '\n']);
                if cfgid[0] == hex2dec('7F') & cfgid[1] == hex2dec('7F'):  # and we have *another* 7F7F
                    valid = 1  # ...ONLY THEN it is valid.
                # end;
            # end;
        except :
            print('Error: seek failed in checkheader')
        # end;
    else:
        fd.seek(-2, cof);
    return valid
# end;

#--------------------------------------
def rd_hdrseg(fd):
    '''Reads a Header
       returns [hdr,nbyte]
    '''
    hdr = Hdr()
    hdr.nbyte = numpy.fromfile(fd, numpy.int16, 1)[0]
    fd.seek(1, cof)
    ndat = numpy.fromfile(fd, numpy.int8, 1)[0]
    hdr.dat_offsets = numpy.fromfile(fd, numpy.int16, ndat)
    nbyte = 4 + ndat * 2;
    return [hdr, nbyte]

#-------------------------------------
def getopt(val, *varargin):
    '''returns one of a list (0=first in varargin, etc.)
    '''
    if val + 1 > len(varargin):
        opt = 'unknown'
    else:
        opt = varargin[val]  # varargin[val + 1];
    # end;
    return opt

def isempty(obj):
    '''
    Checks if the object is empty
    '''
    if isinstance(obj, int):
        if obj == 0:
            return True
        else:
            return False
    elif isinstance(obj, list):
        if len(obj) == 0:
            return True
        else:
            return False

#-------------------------------------
def rd_hdr(fd):
    ''' Read config data
        Changed by Matt Brennan to skip weird stuff at BOF (apparently
        can happen when switching from one flash card to another
        in moored ADCPs).
    return [hdr,pos]
    '''

    cfgid = numpy.fromfile(fd, numpy.uint8, 2)
    nread = 0
    while (cfgid[0] != hex2dec('7F') or cfgid[1] != hex2dec('7F'))  or not checkheader(fd):
        nextbyte = numpy.fromfile(fd, numpy.uint8, 1)[0]
        pos = fd.tell();
        nread = nread + 1
        if isempty(nextbyte):  # End of file
            print('EOF reached before finding valid cfgid')
            hdr = numpy.nan
            return
        # end
        cfgid[2] = cfgid[1]
        cfgid[1] = nextbyte
        if pos % 1000 == 0:
            print('Still looking for valid cfgid at file position %d ...' % pos)
        # end
    # end

    pos = fd.tell() - 2
    if nread > 0:
        print('Junk found at BOF...skipping %d bytes until' % nread)
        print('cfgid=' + dec2hex(cfgid[1]) + dec2hex(cfgid[2]) + ' at file pos %d' % pos)
    # end;

    [hdr, nbyte] = rd_hdrseg(fd)  # ??? was hdr = rd_hdrseg(fd)
    return [hdr, pos]

#
#-------------------------------------
def rd_fixseg(fd):
    '''
        Reads the configuration data from the fixed leader
          return [cfg,nbyte]=
    '''
    # #disp(numpy.fromfile(fd,10,numpy.uint8))
    # #rd.seek(-10,cof);
    cfg = ADCPCfg()
    cfg.name = 'wh-adcp'

    # default - depending on what data blocks are
    # around we can modify this later in rd_buffer.
    cfg.sourceprog = 'instrument';
    cfg.prog_ver = (numpy.fromfile(fd, numpy.uint8, 1) + numpy.fromfile(fd, numpy.uint8, 1) / 100)[0]

    # 8,9,16 - WH navigator
    # 10 -rio grande
    # 15, 17 - NB
    # 19 - REMUS, or customer specific
    # 11- H-ADCP
    # 31 - Streampro
    # 34 - NEMO
    # 50 - WH, no bottom track (built on 16.31)
    # 51 - WH, w/ bottom track
    # 52 - WH, mariner

    if numpy.fix(cfg.prog_ver) == 4 or numpy.fix(cfg.prog_ver) == 5:
        cfg.name = 'bb-adcp'
    elif numpy.fix(cfg.prog_ver) == 8 or numpy.fix(cfg.prog_ver) == 9 or numpy.fix(cfg.prog_ver) == 10 \
            or numpy.fix(cfg.prog_ver) == 16 or numpy.fix(cfg.prog_ver) == 50 or numpy.fix(cfg.prog_ver) == 51  \
            or numpy.fix(cfg.prog_ver) == 52:
        cfg.name = 'wh-adcp'
    elif numpy.fix(cfg.prog_ver) == 14 or numpy.fix(cfg.prog_ver) == 23:  # phase 1 and phase 2
        cfg.name = 'os-adcp'
    else :
        cfg.name = 'unrecognized firmware version'
    # end;

    config = numpy.fromfile(fd, numpy.uint8, 2)  # Coded stuff
    cfg.config = dec2base(config[1], 2, 8) + '-' + dec2base(config[0], 2, 8)
    cfg.beam_angle = getopt((config[1] & 3), 15, 20, 30)
    cfg.numbeams = getopt((config[1] & 16) == 16, 4, 5)
    cfg.beam_freq = getopt((config[0] & 7), 75, 150, 300, 600, 1200, 2400, 38)
    cfg.beam_pattern = getopt((config[0] & 8) == 8, 'concave', 'convex')  # 1=convex,0=concave
    cfg.orientation = getopt((config[0] & 128) == 128, 'down', 'up')  # 1=up,0=down
    cfg.simflag = getopt(numpy.fromfile(fd, numpy.uint8, 1)[0], 'real', 'simulated')  # Flag for simulated data

    fd.seek(1, cof)
    cfg.n_beams = numpy.fromfile(fd, numpy.uint8, 1)[0]
    cfg.n_cells = numpy.fromfile(fd, numpy.uint8, 1)[0]
    cfg.pings_per_ensemble = numpy.fromfile(fd, numpy.uint16, 1)[0]
    cfg.cell_size = numpy.fromfile(fd, numpy.uint16, 1)[0] * .01  # meters
    cfg.blank = numpy.fromfile(fd, numpy.uint16, 1)[0] * .01  # meters
    cfg.prof_mode = numpy.fromfile(fd, numpy.uint8, 1)[0]  #
    cfg.corr_threshold = numpy.fromfile(fd, numpy.uint8, 1)[0]
    cfg.n_codereps = numpy.fromfile(fd, numpy.uint8, 1)[0]
    cfg.min_pgood = numpy.fromfile(fd, numpy.uint8, 1)[0]
    cfg.evel_threshold = numpy.fromfile(fd, numpy.uint16, 1)[0]


    # cfg.time_between_ping_groups=sum(numpy.fromfile(fd,3,numpy.uint8).*[60 1 .01]') # seconds
    ar_mss = numpy.array([60, 1, 0.01]) * numpy.fromfile(fd, numpy.uint8, 3)
    cfg.time_between_ping_groups = ar_mss.sum(axis = 0)  # seconds
    coord_sys = numpy.fromfile(fd, numpy.uint8, 1)[0]  # Lots of bit-mapped info
    cfg.coord = dec2base(coord_sys, 2, 8)
    cfg.coord_sys = getopt(((coord_sys >> 3) & 3), 'beam', 'instrument', 'ship', 'earth')
    cfg.use_pitchroll = getopt((coord_sys & 4) == 4, 'no', 'yes')
    cfg.use_3beam = getopt((coord_sys & 2) == 2, 'no', 'yes')
    cfg.bin_mapping = getopt((coord_sys & 1) == 1, 'no', 'yes')
    cfg.xducer_misalign = numpy.fromfile(fd , numpy.int16, 1)[0] * .01  # degrees
    cfg.magnetic_var = numpy.fromfile(fd, numpy.int16, 1)[0] * .01  # degrees
    cfg.sensors_src = dec2base(numpy.fromfile(fd, numpy.uint8, 1)[0], 2, 8)
    cfg.sensors_avail = dec2base(numpy.fromfile(fd , numpy.uint8, 1)[0], 2, 8)
    cfg.bin1_dist = numpy.fromfile(fd, numpy.uint16, 1)[0] * .01  # meters
    cfg.xmit_pulse = numpy.fromfile(fd, numpy.uint16, 1)[0] * .01  # meters
    cfg.water_ref_cells = numpy.fromfile(fd, numpy.uint8, 2)
    cfg.fls_target_threshold = numpy.fromfile(fd , numpy.uint8, 1)[0]
    fd.seek(1, cof)
    cfg.xmit_lag = numpy.fromfile(fd, numpy.uint16, 1)[0] * .01  # meters
    nbyte = 40


    if numpy.fix(cfg.prog_ver) == 8 or numpy.fix(cfg.prog_ver) == 10 or numpy.fix(cfg.prog_ver) == 16 \
         or numpy.fix(cfg.prog_ver) == 50 or numpy.fix(cfg.prog_ver) == 51 or numpy.fix(cfg.prog_ver) == 52:
        if cfg.prog_ver >= 8.14:  # Added CPU serial number with v8.14
            cfg.serialnum = numpy.fromfile(fd , numpy.uint8, 8)
            nbyte = nbyte + 8;
        # end;


        if cfg.prog_ver >= 8.24:  # Added 2 more bytes with v8.24 firmware
            cfg.sysbandwidth = numpy.fromfile(fd, numpy.uint8, 2)
            nbyte = nbyte + 2
        # end;

        if cfg.prog_ver >= 16.05:  # Added 1 more bytes with v16.05 firmware
            cfg.syspower = numpy.fromfile(fd, numpy.uint8, 1)[0]
            nbyte = nbyte + 1
        # end;



    if cfg.prog_ver >= 16.27:  # Added bytes for REMUS, navigators, and HADCP
        cfg.navigator_basefreqindex = numpy.fromfile(fd, numpy.uint8, 1)[0]
        nbyte = nbyte + 1
        cfg.remus_serialnum = numpy.fromfile(fd, numpy.uint8, 4)
        nbyte = nbyte + 4
        cfg.h_adcp_beam_angle = numpy.fromfile(fd, numpy.uint8, 1)[0]
        nbyte = nbyte + 1
    # end;

    elif numpy.fix(cfg.prog_ver) == 9:
        if cfg.prog_ver >= 9.10:  # Added CPU serial number with v8.14
            cfg.serialnum = numpy.fromfile(fd, numpy.uint8, 8)
            nbyte = nbyte + 8
            cfg.sysbandwidth = numpy.fromfile(fd, numpy.uint8, 2)
            nbyte = nbyte + 2
        # end

    elif numpy.fix(cfg.prog_ver) == 14 or numpy.fix(cfg.prog_ver) == 23:
        cfg.serialnum = numpy.fromfile(fd, numpy.uint8, 8)  # 8 bytes 'reserved'
        nbyte = nbyte + 8

    # end

    # It is useful to have this precomputed.
    transp = numpy.arange(0, cfg.n_cells).transpose().reshape(-1, 1)  # transpose
    cfg.ranges = cfg.bin1_dist + transp[:, 0] * cfg.cell_size
    if cfg.orientation == 'down':  # was ==1
        cfg.ranges = -cfg.ranges
    # end
    return [cfg, nbyte]






#-------------------------------------
def rd_fix(fd):
    '''Read config data'''

    cfgid = numpy.fromfile(fd, numpy.uint16, 1)[0]
    if cfgid != hex2dec('0000'):
        print('Fixed header ID %d incorrect - data corrupted or not a BB/WH raw file?' % cfgid)
    # end;
    [cfg, nbytes] = rd_fixseg(fd)
    return cfg


def get_hdr(fd):
    global cfg

    [hdr, pos] = rd_hdr(fd)
    if not isinstance(hdr, Hdr):  # ~isstruct(hdr):
        ens = -1
        cfg = numpy.nan
        return;
    # end if
    cfg = rd_fix(fd)
    fd.seek(pos, bof)
    # clear global ens
    # global ens
    return [cfg, hdr]

#-----------------------------
def rd_buffer(fd, num_av, debug = False):
    '''
        Reads the ADCP file
        returns [ens,hdr,cfg,pos]=
    '''
    # To save it being re-initialized every time.
    # global hdr,ens

    # A fudge to try and read files not handled quite right.
    global FIXOFFSET
    global SOURCE
    global cfg
    pos = 0

    # If num_av<0 we are reading only 1 element and initializing
    if num_av < 0:
        SOURCE = 0

    # This reinitializes to whatever length of ens we want to average.
    FIXOFFSET = 0

    if num_av < 0:  # | isempty(ens):
        [cfg, hdr] = get_hdr(fd);
        n = abs(num_av)
        ens = Ensamble(n, cfg)
        num_av = abs(num_av)
    else:
        ens = Ensamble(num_av, cfg)
    # end if;

    k = -1;
    while k < num_av - 1:
        # This is in case junk appears in the middle of a file.
        num_search = 6000

        id1 = numpy.fromfile(fd, numpy.uint8, 2)
        search_cnt = 0
        try:
            while search_cnt < num_search and ((id1[0] != hex2dec('7F') or id1[1] != hex2dec('7F')) or not checkheader(fd)):
                search_cnt = search_cnt + 1;
                nextbyte = numpy.fromfile(fd, numpy.uint8, 1)[0]
                if isempty(nextbyte):  # End of file
                    print('EOF reached after %d bytes searched for next valid ensemble start' % search_cnt)
                    ens = -1
                    return
                # end;
                id1[2] = id1[1]
                id1[1] = nextbyte;
                # fprintf([dec2hex(id1(1)) '--' dec2hex(id1(2)) '\n']);
            # end;
        except:
            print('Error reached after %d ensambles read file may be corrupted' % k)

            return [ens, hdr, cfg, pos]

        if search_cnt == num_search:
            print('Searched %d entries...Not a workhorse/broadband file or bad data encountered: -> %x' % (search_cnt, id1))
        elif search_cnt > 0:
            print('Searched %d bytes to find next valid ensemble start' % search_cnt)
        # end

        startpos = fd.tell() - 2  # Starting position.

        # Read the # data types.
        [hdr, nbyte] = rd_hdrseg(fd)
        byte_offset = nbyte + 2
        # # fprintf('# data types = %d\n  ',(length(hdr.dat_offsets)));
        # # fprintf('Blocklen = %d\n  ',hdr.nbyte);

        # Read all the data types.
        for n in range(0, len(hdr.dat_offsets)):
            id = dec2base(numpy.fromfile(fd, numpy.uint16, 1)[0], 16, 4);
            # #   fprintf('ID=%s SOURCE=%d\n',id,SOURCE);

            # handle all the various segments of data. Note that since I read the IDs as a two
            # byte number in little-endian order the high and low bytes are exchanged compared to
            # the values given in the manual.
            #
            winrivprob = 0;

            if id == '0000':  # Fixed leader
                [cfg, nbyte] = rd_fixseg(fd)
                nbyte = nbyte + 2
            elif id == '0080':  # Variable Leader
                k = k + 1
                ens.number[0, k] = numpy.fromfile(fd, numpy.uint16, 1)[0]
                ens.rtc[:, k] = numpy.fromfile(fd, numpy.uint8, 7)
                ens.number[0, k] = ens.number[0, k] + 65536 * numpy.fromfile(fd, numpy.uint8, 1)[0]
                ens.BIT[0, k] = numpy.fromfile(fd, numpy.uint16, 1)[0]
                ens.ssp[0, k] = numpy.fromfile(fd, numpy.uint16, 1)[0]
                ens.depth[0, k] = numpy.fromfile(fd, numpy.uint16, 1)[0] * 0.1  # meters
                ens.heading[0, k] = numpy.fromfile(fd, numpy.uint16, 1)[0] * 0.01  # degrees
                ens.pitch[0, k] = numpy.fromfile(fd, numpy.int16, 1)[0] * 0.01  # degrees
                ens.roll[0, k] = numpy.fromfile(fd, numpy.int16, 1)[0] * 0.01  # degrees
                ens.salinity[0, k] = numpy.fromfile(fd, numpy.int16, 1)[0]  # PSU
                ens.temperature[0, k] = numpy.fromfile(fd, numpy.int16, 1)[0] * 0.01  # Deg C
                ens.mpt[0, k] = sum(numpy.fromfile(fd, numpy.uint8, 3) * [60, 1, 0.01])  # seconds
                ens.heading_std[0, k] = numpy.fromfile(fd, numpy.uint8, 1)[0]  # degrees
                ens.pitch_std[0, k] = numpy.fromfile(fd, numpy.uint8, 1)[0] * 0.1  # degrees
                ens.roll_std[0, k] = numpy.fromfile(fd, numpy.uint8, 1)[0] * 0.1  # degrees
                ens.adc[:, k] = numpy.fromfile(fd, numpy.uint8, 8)
                nbyte = 2 + 40

                if cfg.name == 'bb-adcp':
                    if cfg.prog_ver >= 5.55:
                        fd.seek(15, cof)  # 14 zeros and one byte for number WM4 bytes
                        cent = numpy.fromfile(fd, numpy.uint8, 1)[0]  # possibly also for 5.55-5.58 but
                        ens.rtc[:, k] = numpy.fromfile(fd, numpy.uint8, 7)  # I have no data to test.
                        ens.rtc[0, k] = ens.rtc[0, k] + cent * 100
                        nbyte = nbyte + 15 + 8
                    # end
                elif cfg.name == 'wh-adcp':  # for WH versions.
                    ens.error_status_wd[0, k] = numpy.fromfile(fd, numpy.uint32, 1)[0]
                    nbyte = nbyte + 4

                    if numpy.fix(cfg.prog_ver) == 8 or numpy.fix(cfg.prog_ver) == 10 or numpy.fix(cfg.prog_ver) == 16 \
                    or numpy.fix(cfg.prog_ver) == 50 or numpy.fix(cfg.prog_ver) == 51 or numpy.fix(cfg.prog_ver) == 52:

                        if cfg.prog_ver >= 8.13:  # Added pressure sensor stuff in 8.13
                            fd.seek(2, cof)
                            ens.pressure[0, k] = numpy.fromfile(fd, numpy.uint32, 1)[0]
                            ens.pressure_std[0, k] = numpy.fromfile(fd , numpy.uint32, 1)[0]
                            nbyte = nbyte + 10
                        # end;
                        if cfg.prog_ver >= 8.24:  # Spare byte added 8.24
                            fd.seek(1, cof)
                            nbyte = nbyte + 1
                        # end;

                        if (cfg.prog_ver >= 10.01 and cfg.prog_ver <= 10.99) or cfg.prog_ver >= 16.05:  # Added more fields with century in clock 16.05
                            cent = numpy.fromfile(fd, numpy.uint8, 1)[0]
                            ens.rtc[:, k] = numpy.fromfile(fd, numpy.uint8, 7)
                            ens.rtc[0, k] = ens.rtc[0, k] + int(cent * 100)
                            nbyte = nbyte + 8
                        # end
                    elif numpy.fix(cfg.prog_ver) == 9:
                        fd.seek(2, cof)
                        ens.pressure[0, k] = numpy.fromfile(fd, numpy.uint32, 1)[0]
                        ens.pressure_std[0, k] = numpy.fromfile(fd, numpy.uint32, 1)[0]
                        nbyte = nbyte + 10

                        if cfg.prog_ver >= 9.10:  # Spare byte added 8.24
                            fd.seek(1, cof)
                            nbyte = nbyte + 1;
                        # end;
                    # end;
                elif cfg.name == 'os-adcp':
                    fd.seek(16, cof)  # 30 bytes all set to zero, 14 read above
                    nbyte = nbyte + 16

                    if cfg.prog_ver > 23:
                        fd.seek(2, cof)
                        nbyte = nbyte + 2
                    # end;
                # end;

            elif id == '0100':  # Velocities

                shape = (cfg.n_cells, 4)
                sz = 4 * cfg.n_cells
                # transpose so that vels is a row vector
                ivels = numpy.fromfile(file = fd, dtype = numpy.int16, count = sz).transpose()
                vels = ivels.reshape(shape) * 0.001  # m/s
                ens.east_vel[:, k] = vels[:, 0]
                ens.north_vel[:, k] = vels[:, 1]
                ens.vert_vel[:, k] = vels[:, 2]
                ens.error_vel[:, k] = vels[:, 3]
                nbyte = 2 + 4 * cfg.n_cells * 2

            elif id == '0200':  # Correlations
                shape = (cfg.n_cells, 4)
                sz = 4 * cfg.n_cells
                crr = numpy.fromfile(file = fd, dtype = numpy.uint8, count = sz).reshape(shape)
                ens.corr[:, :, k] = crr
                nbyte = 2 + 4 * cfg.n_cells;

            elif id == '0300':  # Echo Intensities
                shape = (cfg.n_cells, 4)
                sz = 4 * cfg.n_cells
                ints = numpy.fromfile(file = fd, dtype = numpy.uint8, count = sz).reshape(shape)  # .transpose()
                ens.intens[:, :, k] = ints
                nbyte = 2 + 4 * cfg.n_cells

            elif id == '0400':  # Percent good
                shape = (cfg.n_cells, 4)
                sz = 4 * cfg.n_cells
                pct = numpy.fromfile(file = fd, dtype = numpy.uint8, count = sz).reshape(shape)  # .transpose()
                ens.percent[:, :, k] = pct
                nbyte = 2 + 4 * cfg.n_cells

            elif id == '0500':  # Status
                if cfg.name == 'os-adcp':
                    fd.seek(00, cof)
                    nbyte = 2 + 00
                else:
                    # Note in one case with a 4.25 firmware SC-BB, it seems like
                    # this block was actually two bytes short!
                     shape = (cfg.n_cells, 4)
                     sz = 4 * cfg.n_cells
                     st = numpy.fromfile(file = fd, dtype = numpy.uint8, count = sz).reshape(shape)  # .transpose()
                     ens.status[ :, :, k] = st
                     nbyte = 2 + 4 * cfg.n_cells
                # end;

            elif id == '0600':  # Bottom track
                # In WINRIVER GPS data is tucked into here in odd ways, as long
                # as GPS is enabled.

                if SOURCE == 2:
                    fd.seek(2, cof);
                    long1 = numpy.fromfile(fd, numpy.uint16, 1)[0]
                    fd.seek(6 * cof)
                    cfac = 180 / 2 ** 31
                    ens.slatitude[k] = numpy.fromfile(fd, numpy.int32, 1)[0] * cfac
                    if ens.slatitude[0, k] == 0:
                        ens.slatitude[0, k] = numpy.nan
                    # end;
                    # #fprintf('\n k % 8.3f',ens.slatitude(k));
                else:
                    fd.seek(14, cof)  # Skip over a bunch of stuff
                # end;

                ens.bt_range[:, k] = numpy.fromfile(fd, numpy.uint16, 4) * 0.01  #
                ens.bt_vel[:, k] = numpy.fromfile(fd, numpy.int16, 4)
                ens.bt_corr[:, k] = numpy.fromfile(fd, numpy.uint8, 4)  # felipe pimenta aug. 2006
                ens.bt_ampl[:, k] = numpy.fromfile(fd, numpy.uint8, 4)  # "
                ens.bt_perc_good[:, k] = numpy.fromfile(fd, numpy.uint8, 4)  # "
                if SOURCE == 2:
                    rd.seek(2, cof)
                    ens.slongitude[0, k] = (long1 + 65536 * numpy.fromfile(fd, numpy.uint16, 1)[0]) * cfac
                    # #fprintf('\n k % d % 8.3f % f ',long1,ens.slongitude(k),(ens.slongitude(k)/cfac-long1)/65536);
                    if ens.slongitude[0, k] > 180:
                        ens.slongitude[0, k] = ens.slongitude[0, k] - 360
                    # end;
                    if ens.slongitude[0, k] == 0:
                        ens.slongitude[0, k] = numpy.nan
                    # end;
                    rd.seek(16, cof)
                    qual = numpy.fromfile(fd, numpy.uint8, 1)[0]
                    if qual == 0:
                        # # fprintf('qual == % d, % f % f',qual,ens.slatitude(k),ens.slongitude(k));
                        ens.slatitude[0, k] = numpy.nan
                        ens.slongitude[0, k] = NaN
                    # end;
                    rd.seek(71 - 45 - 21, cof)
                else:
                    rd.seek(71 - 45, cof)
                # end;

                nbyte = 2 + 68
                if cfg.prog_ver >= 5.3:  # Version 4.05 firmware seems to be missing these last 11 bytes.
                    rd.seek(78 - 71, cof)
                    ens.bt_range[:, k] = ens.bt_range[:, k] + numpy.fromfile(fd, numpy.uint8, 4) * 655.36
                    nbyte = nbyte + 11

                    if cfg.name == 'wh - adcp':
                        if cfg.prog_ver >= 16.20:  # RDI documentation claims these extra bytes were added in v 8.17
                            rd.seek(4, cof)  # but they don't appear in my 8.33 data - conversation with
                            nbyte = nbyte + 4  # Egil suggests they were added in 16.20
                        # end;
                    # end;
                # end;

                # The raw files produced by VMDAS contain a binary navigation data
                # block.

            elif id == '2000':  # Something from VMDAS.
                cfg.sourceprog = 'VMDAS'
                if SOURCE != 1:
                    print('\n***** Apparently a VMDAS file \n\n')
                # end;
                SOURCE = 1
                utim = numpy.fromfile(fd, numpy.uint8, 4)
                mtime = datenum(utim[3] + utim[4] * 256, utim[2], utim[1])
                ens.smtime[0, k] = mtime + numpy.fromfile(fd, numpy.uint32, 1)[0] / 8640000
                rd.seek(4, cof)  # PC clock offset from UTC
                cfac = 180 / 2 ** 31
                ens.slatitude[0, k] = numpy.fromfile(fd, numpy.int32, [0]) * cfac
                ens.slongitude[0, k] = numpy.fromfile(fd, numpy.int32, 1)[0] * cfac
                ens.emtime[0, k] = mtime + numpy.fromfile(fd, numpy.uint32, 1)[0] / 8640000
                ens.elatitude[0, k] = numpy.fromfile(fd, numpy.int32, 1)[0] * cfac
                ens.elongitude[0, k] = numpy.fromfile(fd, numpy.int32, 1)[0] * cfac
                rd.seek(12, cof)
                ens.flags[0, k] = numpy.fromfile(fd, numpy.uint16, 1)[0]
                rd.seek(6, cof)
                utim = numpy.fromfile(fd, numpy.uint8, 4)
                mtime = datenum(utim[1] + utim[2] * 256, utim[4], utim[3])
                ens.nmtime[0, k] = mtime + numpy.fromfile(fd, numpy.uint32, 1)[0] / 8640000
                # in here we have 'ADCP clock' (not sure how this
                # differs from RTC (in header) and UTC (earlier in this block).
                rd.seek(16, cof)
                nbyte = 2 + 76

            elif id == '2022':  # New NMEA data block from WInRiverII
                cfg.sourceprog = 'WINRIVER2'
                if SOURCE != 2:
                    print('\n***** Apparently a WINRIVER file - Raw NMEA data handler not yet implemented\n\n')
                # end
                SOURCE = 2

                specID = numpy.fromfile(fd, numpy.uint16, 1)[0]
                msgsiz = numpy.fromfile(fd, numpy.int16, 1)[0]
                deltaT = numpy.fromfile(fd, numpy.uchar, 8)
                nbyte = 2 + 12

                rd.seek(msgsiz, cof)
                nbyte = nbyte + msgsiz

                # #    fprintf(' %d ',specID);
                if specID == 100 or specId == 101 or specId == 102 or specID == 103:
                    pass
                # end


            # The following blocks come from WINRIVER files, they aparently contain
            # the raw NMEA data received from a serial port.
            #
            # Note that for WINRIVER files somewhat decoded data is also available
            # tucked into the bottom track block.
            #
            # I've put these all into their own block because RDI's software apparently completely ignores the
            # stated lengths of these blocks and they very often have to be changed. Rather than relying on the
            # error coding at the end of the main block to do this (and to produce an error message) I will
            # do it here, without an error message to emphasize that I am kludging the WINRIVER blocks only!

            elif id == '2100' or id == '2101' or id == '2102' or id == '2103' or id == '2104':

                winrivprob = 1

                if id == '2100':  # $xxDBT  (Winriver addition) 38
                    cfg.sourceprog = 'WINRIVER'
                    if SOURCE != 2:
                        print('\n***** Apparently a WINRIVER file - Raw NMEA data handler not yet implemented\n\n')
                    # end
                    SOURCE = 2
                    strg = numpy.fromfile(fd, numpy.uchar, 38).transpose()
                    nbyte = 2 + 38

                elif id == '2101':  # $xxGGA  (Winriver addition) 94 in manual but 97 seems to work
                    # Except for a winriver2 file which seems to use 77.
                    cfg.sourceprog = 'WINRIVER'
                    if SOURCE != 2:
                        print('\n ***** Apparently a WINRIVER file - Raw NMEA data handler not yet implemented\n\n')
                    # end;
                    SOURCE = 2;

                    asciicode = numpy.fromfile(fd, numpy.uchar, 97).transpose()
                    # initialize the string array
                    strg = numpy.np.zeros((97))
                    for i in asciicode:
                        strg[i] = chr(asciicode[i])
                    nbyte = 2 + 97;
                    l = strfind(strg, '$GPGGA')

                    if not isempty(l):
                        ens.smtime[0, k] = (scanf.sscanf(strg[l + 7:l + 8], '%d') + (scanf.sscanf(strg[l + 9:l + 10], '%d') + scanf.sscanf(strg[l + 11:l + 12], '%d') / 60) / 60) / 24
                    # end
                    #    disp(['->' setstr(str(1:50)) '<-']);

                elif id == '2102':
                    cfg.sourceprog = 'WINRIVER';
                    if SOURCE != 2:
                        print('\n***** Apparently a WINRIVER file - Raw NMEA data handler not yet implemented\n\n')
                    # end;
                    SOURCE = 2;
                    strg = numpy.fromfile(fd, numpy.uchar, 45).transpose()
                    nbyte = 2 + 45
                    # disp(setstr(str));

                elif id == '2103':  # $xxGSA  (Winriver addition) 60
                    cfg.sourceprog = 'WINRIVER';
                    if SOURCE != 2:
                        print('\n ***** Apparently a WINRIVER file - Raw NMEA data handler not yet implemented\n\n')
                    # end
                    SOURCE = 2
                    strg = numpy.fromfile(fd, numpy.uchar, 60).transpose();
                    #  disp(setstr(str));
                    nbyte = 2 + 60

                elif id == '2104':  # xxHDT or HDG (Winriver addition) 38
                    cfg.sourceprog = 'WINRIVER'
                    if SOURCE != 2:
                        print('\n***** Apparently a WINRIVER file - Raw NMEA data handler not yet implemented\n\n')
                    # end
                    SOURCE = 2
                    strg = numpy.fromfile(fd, numpy.uchar, 38).transpose()
                    #      disp(setstr(str));
                    nbyte = 2 + 38
                # end

            elif id == '0701':  # Number of good pings
                rd.seek(4 * cfg.n_cells, cof)
                nbyte = 2 + 4 * cfg.n_cells

            elif id == '0702':  # Sum of squared velocities
                rd.seek(4 * cfg.n_cells, cof)
                nbyte = 2 + 4 * cfg.n_cells

            elif id == '0703':  # Sum of velocities
                rd.seek(4 * cfg.n_cells, cof)
                nbyte = 2 + 4 * cfg.n_cells

                # These blocks were implemented for 5-beam systems

            elif id == '0A00':  # Beam 5 velocity (not implemented)
                rd.seek(cfg.n_cells, cof)
                nbyte = 2 + cfg.n_cells

            elif id == '0301':  # Beam 5 Number of good pings (not implemented)
                rd.seek(cfg.n_cells, cof)
                nbyte = 2 + cfg.n_cells

            elif id == '0302':  # Beam 5 Sum of squared velocities (not implemented)
                rd.seek(cfg.n_cells, cof)
                nbyte = 2 + cfg.n_cells

            elif id == '0303':  # Beam 5 Sum of velocities (not implemented)
                rd.seek(cfg.n_cells, cof)
                nbyte = 2 + cfg.n_cells

            elif id == '020C':  # Ambient sound profile (not implemented)
                rd.seek(4, cof)
                nbyte = 2 + 4

            elif id == '3000':  # Fixed attitude data format for OS-ADCPs (not implemented)
                rd.seek(32, cof)
                nbyte = 2 + 32

            else:
                # This is pretty idiotic - for OS-ADCPs (phase 2) they suddenly decided to code
                # the number of bytes into the header ID word. And then they don't really
                # document what they did! So, this is cruft of a high order, and although
                # it works on the one example I have - caveat emptor....
                #
                # Anyway, there appear to be codes 0340-03FC to deal with. I am not going to
                # decode them but I am going to try to figure out how many bytes to
                # skip.
                if id[1:2] == '30':
                    # I want to count the number of 1s in the middle 4 bits of the
                    # 2nd two bytes.
                    nflds = sum(dec2base((hex2dec(id[3:4]) & hex2dec('3C')), 2) == '1')
                    # I want to count the number of 1s in the highest 2 bits of byte 3
                    dfac = sum(dec2base((hex2dec(id[3]) & hex2dec('C')), 2) == '1')
                    rd.seek(12 * nflds * dfac, cof)
                    nbyte = 2 + 12 * nflds * dfac

                else:
                    print('Unrecognized ID code: %s\n' % id)
                    nbyte = 2
                # end;
                # # ens=-1;
                # # return;
            # end; BIG IF

            # here I adjust the number of bytes so I am sure to begin
            # reading at the next valid offset. If everything is working right I shouldn't have
            # to do this but every so often firware changes result in some differences.

            # #fprintf('#bytes is %d, original offset is %d\n',nbyte,byte_offset);
            byte_offset = byte_offset + nbyte;

            if n < len(hdr.dat_offsets) - 1:
                if hdr.dat_offsets[n + 1] != byte_offset:
                    if not winrivprob and debug:
                        print('%s: Adjust location by %d\n' % (id, hdr.dat_offsets[n + 1] - byte_offset))
                    # end;
                    fd.seek(hdr.dat_offsets[n + 1] - byte_offset, cof);
                # end;
                byte_offset = hdr.dat_offsets[n + 1]
            else:
                if hdr.nbyte - 2 != byte_offset:
                    if not winrivprob and debug:
                        print('%s: Adjust location by %d\n' % (id, hdr.nbyte - 2 - byte_offset))
                    # end;
                    fd.seek(hdr.nbyte - 2 - byte_offset, cof)
                # end
                byte_offset = hdr.nbyte - 2
            # end
        # end; BIG for loop

        # Now at the end of the record we have two reserved bytes, followed
        # by a two-byte checksum = 4 bytes to skip over.

        readbytes = fd.tell() - startpos
        offset = (hdr.nbyte + 2) - byte_offset  # The 2 is for the checksum

        if offset != 4 and FIXOFFSET == 0:
            print('\n*****************************************************\n')
            if fd.feof():
                print(' EOF reached unexpectedly - discarding this last ensemble\n')
                ens = -1
            else:
                print('Adjust location by %d (readbytes=%d, hdr.nbyte=%d)\n', offset, readbytes, hdr.nbyte)
                print(' NOTE - If this appears at the beginning of the read, it is\n')
                print('        is a program problem, possibly fixed by a fudge\n')
                print('        PLEASE REPORT TO rich@eos.ubc.ca WITH DETAILS!!\n\n')
                print('      -If this appears at the end of the file it means\n')
                print('       The file is corrupted and only a partial record has  \n')
                print('       has been read\n')
            # end;
            print('******************************************************\n')
            FIXOFFSET = offset - 4
        # end;
        fd.seek(4 + FIXOFFSET, cof)

        # An early version of WAVESMON and PARSE contained a bug which stuck an additional two
        # bytes in these files, but they really shouldn't be there
        # if cfg.prog_ver>=16.05,
        #      rd.seek(2,cof);
        # end;

    # end;

    # Blank out stuff bigger than error velocity
    # big_err=abs(ens.error_vel)>.2;
    big_err = 0;

    # Blank out invalid data
    for i in range(0, len(ens.east_vel[:, 0])) :  # length of first column
        for j in range(0, len(ens.east_vel[0, :])):
            if (ens.east_vel[i, j] == -32.768 or big_err):
                ens.east_vel[i, j] = numpy.nan
            if (ens.north_vel[i, j] == -32.768 or big_err):
                ens.north_vel[i, j] = numpy.nan
            if (ens.vert_vel[i, j] == -32.768 or big_err):
                ens.vert_vel[i, j] = numpy.nan
            if (ens.error_vel[i, j] == -32.768 or big_err):
                ens.error_vel[i, j] = numpy.nan

    return [ens, hdr, cfg, pos]


def  readRawBinADCP(name, *varargin):
    '''
        See descritption at the top of the file
    '''
    num_av = 5  # Block filtering and decimation parameter (# ensembles to block together).
    nens = -1  # Read all ensembles.
    century = 2000  # ADCP clock does not have century prior to firmware 16.05.
    vels = 'no'  # Default to simple averaging
    printinfo = False
    debug = False

    aidx = 0
    for arg in varargin:
        print(" arg[%d]: %s" % (aidx, arg))
        aidx += 1


    nv = len(varargin)
    lv = 0
    # Read optional args
    while lv < nv:
        if lv == 0 and not isinstance(varargin[lv], str):  # ~isstr(varargin{1}),
            num_av = varargin[lv]  # Block filtering and decimation parameter (# ensembles to block together).
            lv += 1
        elif lv >= 1 and not isinstance(varargin[lv], str):
            nens = varargin[lv]
            lv += 1
        elif lv >= 2:
            varvalue = varargin[lv][0:3]
            if  varvalue == 'bas':
                century = varargin[lv + 1]
                lv += 2
            elif  varvalue == 'des':
                if isinstance(varargin[lv], str):
                    if varargin[lv] == 'no':
                        vels = 'no'
                    else:  # yes
                        vels = [0.3, 0.3, 0.3]
                    # end
                else:
                    vels = varargin[2]
                # end
                lv += 2
            elif varvalue == "inf":
                if isinstance(varargin[lv], str):
                    if varargin[lv + 1] == 'yes':
                        printinfo = True
                    else:
                        printinfo = False
                lv += 2
            elif varvalue == "deb":
                if isinstance(varargin[lv], str):
                    if varargin[lv + 1] == 'yes':
                        debug = True
                    else:
                        debug = False
                lv += 2
            else:
                raise Exception('Unknown command line option ' + varargin[lv])
            # end

    # end for

    # Check file information first
    try:
        # returns a tuple (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime)
        st = os.stat(name)
    except IOError:
        print("failed to get information about", name)
    else:
        print("file size:%d" % st.st_size)
        print("file modified: %s" % time.asctime(time.localtime(st.st_mtime)))


    print('\n Opening file %s \n\n' % name)
    # fd=fopen(name,'r','ieee-le');
    fd = open(name, 'rb')
    # Read first ensemble to initialize parameters

    if printinfo == True:
        [ens, hdr, cfg, pos] = rd_buffer(fd, -2, debug);  # Initialize and read first two records
        if not isinstance(ens, Ensamble) and ens == -1:
            print('No Valid data found!')
            adcp = []
            return
        # end;
        fd.seek(pos, bof);  # Rewind

        if (cfg.prog_ver < 16.05 and cfg.prog_ver > 5.999) or cfg.prog_ver < 5.55:
            print('***** Assuming that the century begins year %d (info not in this firmware version) \n\n' % century)
        else :
            century = 0  # century included in clock
        # end;

        # dats = datenum(century + ens.rtc(1, :), ens.rtc(2, :), ens.rtc(3, :), ens.rtc(4, :), ens.rtc(5, :), ens.rtc(6, :) + ens.rtc(7, :) / 100);
        datetm1 = datetime(century + int(ens.rtc[0, 0]), int(ens.rtc[1, 0]), int(ens.rtc[2, 0]), int(ens.rtc[3, 0]), int(ens.rtc[4, 0]), int(ens.rtc[5, 0] + ens.rtc[6, 0] / 100))
        datetm2 = datetime(century + int(ens.rtc[0, 1]), int(ens.rtc[1, 1]), int(ens.rtc[2, 1]), int(ens.rtc[3, 1]), int(ens.rtc[4, 1]), int(ens.rtc[5, 1] + ens.rtc[6, 1] / 100))
        datenm1 = date2num(datetm1)
        datenm2 = date2num(datetm2)
        datest1 = datetm1.strftime("%Y-%m-%d %H:%M:%S")
        datest2 = datetm2.strftime("%Y-%m-%d %H:%M:%S")

        t_int = abs(datetm1 - datetm2)

        print('Record begins at %s\n' % datest1)
        print('Ping interval appears to be:\n \tdays:%d\n\tseconds:%d\n\tmicroseconds:%d\n\n' % (t_int.days, t_int.seconds, t_int.microseconds))
    else:
        [cfg, hdr] = get_hdr(fd);

    # Estimate number of records (since I don't feel like handling EOFs correctly,
    # we just don't read that far!)


    # Now, this is a puzzle - it appears that this is not necessary in
    # a firmware v16.12 sent to me, and I can't find any example for
    # which it *is* necessary so I'm not sure why its there. It could be
    # a leftoever from dealing with the bad WAVESMON/PARSE problem (now
    # fixed) that inserted extra bytes.
    # ...So its out for now.
    # if cfg.prog_ver>=16.05, extrabytes=2; else extrabytes=0; end; # Extra bytes
    extrabytes = 0

    nensinfile = numpy.fix(st.st_size / (hdr.nbyte + 2 + extrabytes));
    print('\nEstimating %d ensembles in this file\n' % nensinfile)

    if len(nens) == 1:
        if nens == -1:
            nens = nensinfile
        # end;
        print('   Reading %d ensembles, reducing by a factor of %d\n' % (nens[0], num_av))
    else :
        print('   Reading ensembles %d-%d, reducing by a factor of %d\n' % (nens[0], nens[1], num_av))
        skipb = (hdr.nbyte + 2 + extrabytes) * (nens[0] - 1)
        fd.seek(skipb, cof)
        nens = numpy.diff(nens, n = 1, axis = 0) + 1
    # end;

    # Number of records after averaging.

    n = int(numpy.fix(nens / num_av));
    print('Final result %d values\n' % n)

    if num_av > 1:
        if isinstance(vels, str):
            print('\n Simple mean used for ensemble averaging\n')
        else:
            print('\n Averaging after outlier rejection with parameters [%f %f %f]\n', vels)
        # end
    # end


    # Structure to hold all ADCP data
    # Note that I am not storing all the data contained in the raw binary file, merely
    # things I think are useful.
    adcp = ADCPData(n, cfg)

    # BH added this for persistance
    adcp.num_rec = n

    # Calibration factors for backscatter data

    # clear global ens ??

    # Loop for all records
    for k in range(0, n - 1):

      # Gives display so you know something is going on...
        if printinfo == True:
            if k % 50 == 0:
                print('\n%d' % (k * num_av), end=' ')
            # end;
            print('.', end=' ')

      # Read an ensemble

        [ens, hdr, cfg, pos] = rd_buffer(fd, num_av, debug)

        if not isinstance(ens, Ensamble):  # If aborting...
            print('Only %d records found..suggest re-running RDRADCP using this parameter\n' % (k - 1) * num_av)
            # print '(If this message preceded by a POSSIBLE PROGRAM PROBLEM message, re-run using %d)\n' % (k - 1) * num_av - 1
            break
        # end;

        # dats = datenum(century + ens.rtc(1, :), ens.rtc(2, :), ens.rtc(3, :), ens.rtc(4, :), ens.rtc(5, :), ens.rtc(6, :) + ens.rtc(7, :) / 100);
        datetm = []
        datenm = []
        for ii in range(0, ens.rtc[:, 0].ndim):
            datetm.insert(ii, datetime(int(century + ens.rtc[0, ii]), int(ens.rtc[1, ii]), int(ens.rtc[2, ii]), int(ens.rtc[3, ii]), int(ens.rtc[4, ii]), int(ens.rtc[5, 0] + ens.rtc[6, ii] / 100)))
            datenm.insert(ii, date2num(datetm[ii]))
        # end
        adcp.mtime[0, k] = numpy.median(datenm)
        adcp.number[0, k] = ens.number[0, 0]
        adcp.heading[0, k] = numpy.mean(ens.heading)
        adcp.pitch[0, k] = numpy.mean(ens.pitch)
        adcp.roll[0, k] = numpy.mean(ens.roll)
        adcp.heading_std[0, k] = numpy.mean(ens.heading_std)
        adcp.pitch_std[0, k] = numpy.mean(ens.pitch_std)
        adcp.roll_std[0, k] = numpy.mean(ens.roll_std)
        adcp.depth[0, k] = numpy.mean(ens.depth)
        adcp.temperature[0, k] = numpy.mean(ens.temperature)
        adcp.salinity[0, k] = numpy.mean(ens.salinity)
        adcp.pressure[0, k] = numpy.mean(ens.pressure)
        adcp.pressure_std[0, k] = numpy.mean(ens.pressure_std)

        if isinstance(vels, str):
            adcp.east_vel[:, k] = nmean(ens.east_vel , 2)[:, 0]
            adcp.north_vel[:, k] = nmean(ens.north_vel, 2)[:, 0]
            adcp.vert_vel[:, k] = nmean(ens.vert_vel , 2)[:, 0]
            adcp.error_vel[:, k] = nmean(ens.error_vel, 2)[:, 0]
        else:
            adcp.east_vel[:, k] = nmedian(ens.east_vel  , vels[0], 2)[:, 0]
            adcp.north_vel[:, k] = nmedian(ens.north_vel, vels[0], 2)[:, 0]
            adcp.vert_vel[:, k] = nmedian(ens.vert_vel  , vels[1], 2)[:, 0]
            adcp.error_vel[:, k] = nmedian(ens.error_vel, vels[2], 2)[:, 0]
        # end

        adcp.corr[:, :, k] = nmean(ens.corr, 2)  # [:, 0]#nmean(ens.corr, 2)        # added correlation RKD 9/00
        adcp.status[:, :, k] = nmean(ens.status, 2)  # [:, 0]#nmean(ens.status[0], 3)
        # print nmean(ens.intens, 2)
        # print ens.intens
        # print numpy.mean(ens.intens, 2)
        adcp.intens[:, :, k] = nmean(ens.intens, 2)  # nmean(ens.intens[0], 3)
        adcp.perc_good[:, :, k] = nmean(ens.percent, 2)  # [:, 0]#nmean(ens.percent[0], 3)  # felipe pimenta aug. 2006 - BAD BAD code


        adcp.bt_range[:, k] = nmean(ens.bt_range, 1)  # nmean(ens.bt_range, 2)
        adcp.bt_vel[:, k] = nmean(ens.bt_vel, 1)  # nmean(ens.bt_vel, 2)

        adcp.bt_corr[:, k] = nmean(ens.bt_corr, 1)  # nmean(ens.bt_corr, 2)          # felipe pimenta aug. 2006
        adcp.bt_ampl[:, k] = nmean(ens.bt_ampl, 1)  # nmean(ens.bt_ampl, 2)          #  "
        adcp.bt_perc_good[:, k] = nmean(ens.bt_perc_good, 1)  # nmean(ens.bt_perc_good, 2)#  "

        if cfg.sourceprog == 'WINRIVER':
            adcp.nav_mtime[0, k] = numpy.mean(ens.smtime)  # nmean(ens.smtime)
            adcp.nav_longitude[0, k] = numpy.mean(ens.slongitude)  # nmean(ens.slongitude)
            adcp.nav_latitude[0, k] = numpy.mean(ens.slatitude)  # nmean(ens.slatitude)
        elif cfg.sourceprog == 'VMDAS':
            adcp.nav_smtime[0, k] = ens.smtime[0]
            adcp.nav_emtime[0, k] = ens.emtime[0]
            adcp.nav_slatitude[0, k] = ens.slatitude[0]
            adcp.nav_elatitude[0, k] = ens.elatitude[0]
            adcp.nav_slongitude[0, k] = ens.slongitude[0]
            adcp.nav_elongitude[0, k] = ens.elongitude[0]
            adcp.nav_mtime[0, k] = nmean(ens.nmtime)
        # end;
    # end;

    print('\n')
    print('Read to byte %d in a file of size %d bytes' % (fd.tell(), st.st_size))
    if fd.tell() + hdr.nbyte < st.st_size:
        print('-->There may be another %d ensembles unread' % numpy.fix((st.st_size - fd.tell()) / (hdr.nbyte + 2)))
    # end

    fd.close();
    adcp.goodbins = int(numpy.round((adcp.depth[0][1000] - adcp.config.bin1_dist) / adcp.config.cell_size))
    return [adcp, cfg, ens, hdr]


#--------------------------------------
def nmean(x, dim = -1):
    # R_NMEAN Computes the mean of matrix ignoring NaN
    #         values
    #   R_NMEAN(X,DIM) takes the mean along the dimension DIM of X.
    #

    # kk = numpy.isfinite(x);
    # x[ not kk] = 0
    x[numpy.isnan(x)] = 0  # set to zero NAN values
    x[numpy.isinf(x)] = 0  # set to 0 infinite values
    if dim == -1:
        # Determine which dimension SUM will use, find nonzero elements
        dim = math.min((x.shape != 1).nonzero())
        if isempty(dim):
             dim = 1
        # end
    # end

    # if dim > length(size(x))
    if dim > x.ndim:
        y = x  # For matlab 5.0 only!!! Later versions have a fixed 'sum'
    else:
        kk = numpy.copy(x)
        kk[numpy.isfinite(x)] = 1
        ndat = kk.sum(axis = dim)  # kk.sum(axis = dim - 1)
        # indat = numpy.copy(ndat[:, ndat == 0])
        indat = numpy.where(ndat == 0)[0]
        # ndat[:, ndat == 0] = 1  # If there are no good data then it doesn't matter what
        ndat[indat] = 1  # If there are no good data then it doesn't matter what
                                  # we average by - and this avoid div-by-zero warnings.
        y = x.sum(axis = dim) / ndat  # x.sum(axis = dim - 1) / ndat
        # y[:, indat == 0] = numpy.nan
        y[indat] = numpy.nan
    # end
    return y

#--------------------------------------
def nmedian(x, window, dim = -1):
    '''Copied from median but with handling of NaN different.

    '''

    if dim == -1:
        dim = math.min((x.size != 1).nonzero())
        if isempty(dim):
            dim = 1
        # end
    # end

    siz = x.shape
    n = x.shape[dim - 1]  # size(x, dim);

    # Permute and reshape so that DIM becomes the row dimension of a 2-D array
    # order = [dim:max(length(size(x)), dim) 1:dim - 1]
    # create an array to help permute
    # x = reshape(permute(x, order), n, prod(siz) / n);
    x = x.transpose()
    # Sort along first dimension. Is it column?
    # x = sort(x, 1);
    x = numpy.sort(x, axis = 0)
    (n1, n2) = x.shape
    if n1 == 1:
        y = x
    else:
        if n2 == 1:
            kk = numpy.sum(numpy.isfinite(x), 1)
            if kk > 0:
                x1 = x((kk[0] - 1) / 2 + 1)
                x2 = x((kk[0] / 2) + 1)
                x[numpy.abs(x - (x1 + x2) / 2) > window] = numpy.nan
            # end
            x = numpy.sort(x, axis = 0)
            kk = numpy.sum(numpy.isfinite(x), 0)
            x[numpy.isnan(x)] = 0
            y = numpy.nan
            if kk > 0:
                y = numpy.sum(x, axis = 0) / kk
            # end
        else:
            kk = numpy.sum(numpy.isfinite(x), 0)
            ll = kk < n1 - 2
            kk[ll] = 0;
            x[:, ll] = numpy.nan
            x1 = x[(kk[0] - 1) / 2, ]
            x2 = x[(kk[0] / 2), ]
            x[numpy.abs(x - numpy.ones((n1, 1), int) * (x1 + x2) / 2) > window] = numpy.nan
            x = numpy.sort(x, axis = 0)
            kk = numpy.sum(numpy.isfinite(x), 0)
            x[numpy.isnan(x)] = 0
            y = numpy.nan + numpy.ones((1, n2))
            if numpy.any(kk):
                y[:, kk > 0] = numpy.sum(x[:, kk > 0], axis = 0) / kk[kk > 0]
            # end
        # end
    # end

    # Permute and reshape back
    # siza = numpy.zeros((siz))
    sz = [ (siz[i]) for i in range(0, len(siz))]
    sz[dim - 1] = 1
    siz = ()
    for i in range(0, len(sz)):
        siz = siz + (sz[i],)
    y = y.transpose().reshape(siz)

    return y
















