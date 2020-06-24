import csv
import numpy as np
from rdradcp.readRawBinADCP import ADCPData
from rdradcp.readRawBinADCP import ADCPCfg
from utools.isnumber import isnumber
from matplotlib.dates import date2num, num2date

class RdiDataWriter(object):
    """
    1) Write an adcp object to files and
    2) Read an adcp object from files (must have the same name and reside in the same directory

    adcp =
        self.goodbins = 0
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


        cfg =
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
    """


    #------------------------------------------------------------------------
    def __init__(self, path, adcp=None, n=0, fname=None):

        self.path = path
        if adcp is None:
            cfg = self.readCfg(fname+'.ctl')
            self.adcp = ADCPData(n, cfg)
        else:
            self.adcp = adcp



      #------------------------------------------------------------------------
    def readAdcp(self, fname, delimiter=',', interval=None):
        '''

        :param fname:
        :param delimiter:
        :return:
        '''
        adcp = self.adcp

        self.readData(fname + ".dat", adcp, delimiter)
        adcp.mtime, adcp.east_vel= self.readBins(fname+'.ve', adcp.config.n_cells, adcp.num_rec, delimiter, interval)
        adcp.mtime, adcp.north_vel = self.readBins(fname + '.vn', adcp.config.n_cells, adcp.num_rec, delimiter, interval)
        adcp.mtime, adcp.vert_vel = self.readBins(fname + '.vu', adcp.config.n_cells, adcp.num_rec, delimiter, interval)
        adcp.mtime, adcp.error_vel = self.readBins(fname + '.err', adcp.config.n_cells, adcp.num_rec, delimiter, interval)

        return adcp


        # ------------------------------------------------------------------------

    def readBins(self, fname, n_cells,  num_rec, delimiter=',', interval=None):
        '''

        :param fname:
        :param delimiter:
        :return:
        '''
        ifile = open(self.path + '/' + fname, "rt")
        reader = csv.reader(ifile, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)

        date = np.zeros((1, num_rec))
        data = np.zeros((n_cells, num_rec))

        rowno = 0
        if interval is not None:
            start_num = interval[0]
            end_num = interval[1]

        for row in reader:
            if rowno == 0: #skip header
                rowno += 1
                continue
            #     colno = 0
            #     for col in row: #create the matrix for all  bins
            #         if colno == 0:
            #             colno += 1
            #             continue
            #         data.append([])
            #         colno += 1
            #     rowno += 1
            #     continue #skip the header line

            colno=0
            for col in row:
                if colno == 0: #skip the string date
                    colno += 1
                    continue
                if colno == 1:
                    fcol = float(col)
                    if interval is not None:
                        if start_num > fcol or end_num < fcol:
                            # this row does not qualify
                            break
                    date[0, rowno-1] = fcol
                else:
                    data[colno-2, rowno-1] = col
                colno += 1
            rowno += 1
        return date, data


    #------------------------------------------------------------------------
    def readCfg(self, fname, delimiter=','):
        '''

        :param fname:
        :return:
        '''
        ifile = open(self.path + '/' + fname, 'r', encoding='utf-8', newline='')
        reader = csv.reader(ifile, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_NONE)
        obj = ADCPCfg()

        for line in reader:
            token = line[0]

            if token == 'ranges':
                # attr = np.fromstring(line[1], dtype=float, sep=',')
                rngs = line[1]
                unicode_line = rngs.translate({ord(c): None for c in '[]"'})
                attr = np.fromstring(unicode_line, dtype=float, sep=' ')
            else:
                if isnumber(line[1]):
                    if '.' in line[1] or 'e' in line[1] or 'E' in line[1]:
                        attr = float(line[1])
                    else:
                        attr = int(line[1])
                else:
                    attr = line[1]
            setattr(obj, token, attr)

        return obj



    # ------------------------------------------------------------------------
    @staticmethod
    def write_binfile(writer, date, velocity, delimiter=','):
        idx = 0
        numdat = []
        prev = 0
        prevtxt = ''
        bins = velocity.shape[0]
        newrow = ['date']
        for i in range(0, velocity.shape[0]):
            newrow.append('bin%d' % i)
        writer.writerow(newrow)

        for i in range(0, velocity.shape[1]):
            onedate = date[0][i]
            try:
                dt = num2date(onedate)

                # Mike 3 has ASCII exported in format 2017-05-31 23:00:00
                datestr = dt.strftime("%Y-%m-%d %H:%M:%S")
                newrow = [datestr]
                newrow.append(onedate)

                for j in range(0, velocity.shape[0]):
                    newrow.append(velocity[j, i])

                writer.writerow(newrow)
            except:
                pass #skip bad date values
        print("Done writing binfile")


    # ------------------------------------------------------------------------
    def writeBins(self, fname, date, data, delimiter):
        '''

        :param fname:
        :param date:
        :param data:
        :param delft3d:
        :return:
        '''


        ofile = open(self.path + '/' + fname, "wt")
        writer = csv.writer(ofile, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
        self.write_binfile(writer, date, data, delimiter)
        ofile.close()


    # ------------------------------------------------------------------------
    def readData(self, fname, adcp, delimiter=',', interval=None):
        '''

        :param fname:
        :param adcp:
        :param delimiter:
        :return:
        '''
        ifile = open(self.path + '/' + fname, "rt")
        reader = csv.reader(ifile, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for row in reader:

            if row[0] == "goodbins":
                adcp.goodbins = int(row[1])

            if row[0] == "num_rec":
                adcp.num_rec = int(row[1])

        ifile.close()


    # ------------------------------------------------------------------------
    def writeData(self, fname, delimiter=','):
        '''

        :param fname:
        :param cfg:
        :param delimiter:
        :return:
        '''
        ofile = open(self.path + '/' + fname, "wt")
        writer = csv.writer(ofile, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)

        if self.adcp.goodbins != 0:
            writer.writerow(["goodbins", int(self.adcp.goodbins)])

        writer.writerow(["num_rec", int(self.adcp.num_rec)])

        ofile.close()


    #------------------------------------------------------------------------
    def writeCfg(self, fname, cfg, delimiter=','):
        '''

        :param fname:
        :param cfg:
        :param delimiter:
        :return:
        '''
        ofile = open(self.path + '/' + fname, "wt")
        writer = csv.writer(ofile, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)

        if cfg.prog_ver != '':
            writer.writerow(["prog_ver", cfg.prog_ver])
        else:
            writer.writerow(["prog_ver", ""])

        if cfg.config is None:
            writer.writerow(["config", "None"])
        else:
            writer.writerow(["config", cfg.config])

        if cfg.beam_angle is None:
            writer.writerow(["beam_angle", "None"])
        else:
            writer.writerow(["beam_angle", cfg.beam_angle])

        if cfg.numbeams is None:
            writer.writerow(["numbeams", "None"])
        else:
            writer.writerow(["numbeams", cfg.numbeams])

        if cfg.beam_freq is None:
            writer.writerow(["beam_freq", "None"])
        else:
            writer.writerow(["beam_freq", cfg.beam_freq])

        if cfg.beam_pattern is None:
            writer.writerow(["beam_pattern", "None"])
        else:
            writer.writerow(["beam_pattern", cfg.beam_pattern])

        if cfg.orientation is None:
            writer.writerow(["orientation", "None"])
        else:
            writer.writerow(["orientation", cfg.orientation])

        if cfg.simflag is None:
            writer.writerow(["simflag", "None"])
        else:
            writer.writerow(["simflag", cfg.simflag])

        if cfg.n_beams is None:
            writer.writerow(["n_beams", "None"])
        else:
            writer.writerow(["n_beams", cfg.n_beams])

        if cfg.n_cells is None:
            writer.writerow(["n_cells", "None"])
        else:
            writer.writerow(["n_cells", cfg.n_cells])

        if cfg.pings_per_ensemble is None:
            writer.writerow(["pings_per_ensemble", "None"])
        else:
            writer.writerow(["pings_per_ensemble", cfg.pings_per_ensemble])

        if cfg.cell_size is None:
            writer.writerow(["cell_size", "None"])
        else:
            writer.writerow(["cell_size", cfg.cell_size])

        if cfg.blank is None:
            writer.writerow(["blank", "None"])
        else:
            writer.writerow(["blank", cfg.blank])

        if cfg.prof_mode is None:
            writer.writerow(["prof_mode", "None"])
        else:
            writer.writerow(["prof_mode", cfg.prof_mode])

        if cfg.corr_threshold is None:
            writer.writerow(["corr_threshold", "None"])
        else:
            writer.writerow(["corr_threshold", cfg.corr_threshold])

        if cfg.n_codereps is None:
            writer.writerow(["n_codereps", "None"])
        else:
            writer.writerow(["n_codereps", cfg.n_codereps])

        if cfg.min_pgood is None:
            writer.writerow(["min_pgood", "None"])
        else:
            writer.writerow(["min_pgood", cfg.min_pgood])

        if cfg.evel_threshold is None:
            writer.writerow(["evel_threshold", "None"])
        else:
            writer.writerow(["evel_threshold", cfg.evel_threshold])

        if cfg.time_between_ping_groups is None:
            writer.writerow(["time_between_ping_groups", "None"])
        else:
            writer.writerow(["time_between_ping_groups", cfg.time_between_ping_groups])

        if cfg.coord is None:
            writer.writerow(["coord", "None"])
        else:
            writer.writerow(["coord", cfg.coord])

        if cfg.coord_sys is None:
            writer.writerow(["coord_sys", "None"])
        else:
            writer.writerow(["coord_sys", cfg.coord_sys])

        if cfg.use_pitchroll is None:
            writer.writerow(["use_pitchroll", "None"])
        else:
            writer.writerow(["use_pithcroll", cfg.use_pitchroll])

        if cfg.use_3beam is None:
            writer.writerow(["use_3beam", "None"])
        else:
            writer.writerow(["use_3beam", cfg.use_3beam])

        if cfg.bin_mapping is None:
            writer.writerow(["bin_mapping", "None"])
        else:
            writer.writerow(["bin_mapping", cfg.bin_mapping])

        if cfg.xducer_misalign is None:
            writer.writerow(["xducer_misalign", "None"])
        else:
            writer.writerow(["xducer_misalign", cfg.xducer_misalign])

        if cfg.magnetic_var is None:
            writer.writerow(["magnetic_var", "None"])
        else:
            writer.writerow(["magnetic_var", cfg.magnetic_var])

        if cfg.sensors_src is None:
            writer.writerow(["sensors_src", "None"])
        else:
            writer.writerow(["sensors_src", cfg.sensors_src])

        if cfg.sensors_avail is None:
            writer.writerow(["sensors_avail", "None"])
        else:
            writer.writerow(["sensors_avail", cfg.sensors_avail])

        if cfg.bin1_dist is None:
            writer.writerow(["bin1_dist", "None"])
        else:
            writer.writerow(["bin1_dist", cfg.bin1_dist])

        if cfg.xmit_pulse is None:
            writer.writerow(["xmit_pulse", "None"])
        else:
            writer.writerow(["xmit_pulse", cfg.xmit_pulse])

        if cfg.water_ref_cells is None:
            writer.writerow(["water_ref_cells", "None"])
        else:
            writer.writerow(["water_ref_cells", cfg.water_ref_cells])

        if cfg.fls_target_threshold is None:
            writer.writerow(["fls_target_threshold", "None"])
        else:
            writer.writerow(["fls_target_threshold", cfg.fls_target_threshold])

        if cfg.xmit_lag is None:
            writer.writerow(["xmit_lag", "None"])
        else:
            writer.writerow(["xmit_lag", cfg.xmit_lag])

        if cfg.serialnum is None:
            writer.writerow(["serialnum", "None"])
        else:
            writer.writerow(["serialnum", cfg.serialnum])

        if cfg.sysbandwidth is None:
            writer.writerow(["sysbandwidth", "None"])
        else:
            writer.writerow(["sysbandwidth", cfg.sysbandwidth])

        if cfg.syspower is None:
            writer.writerow(["syspower", "None"])
        else:
            writer.writerow(["syspower", cfg.syspower])

        if cfg.navigator_basefreqindex is None:
            writer.writerow(["navigator_basefreqindex", "None"])
        else:
            writer.writerow(["navigator_basefreqindex", cfg.navigator_basefreqindex])

        if cfg.remus_serialnum is None:
            writer.writerow(["remus_serialnum", "None"])
        else:
            writer.writerow(["remus_serialnum", cfg.remus_serialnum])

        if cfg.h_adcp_beam_angle is None:
            writer.writerow(["h_adcp_beam_angle", "None"])
        else:
            writer.writerow(["h_adcp_beam_angle", cfg.h_adcp_beam_angle])

        if cfg.ranges is None:
            writer.writerow(["ranges", "None"])
        else:
            writer.writerow(["ranges", cfg.ranges])

        ofile.close()
        print("Done writing " + fname + " file")


    #------------------------------------------------------------------------
    def writeAdcp(self, fname, adcp=None, delimiter=','):
        '''

        :param fname:
        :param date:
        :param adcp:
        :param delimiter:
        :return:
        '''
        if adcp is None:
            adcp = self.adcp
        #self.writeCfg(fname+".ctl", adcp.config, delimiter=delimiter)
        for index, item in enumerate(adcp):
            if type(item).__name__ == 'ADCPCfg':
                self.writeCfg(fname+".ctl", item, delimiter=delimiter)
        self.writeBins(fname+".ve", adcp.mtime, adcp.east_vel, delimiter)
        self.writeBins(fname + ".vn", adcp.mtime, adcp.north_vel, delimiter)
        self.writeBins(fname + ".vu",  adcp.mtime, adcp.vert_vel,  delimiter)
        self.writeBins(fname + ".err", adcp.mtime, adcp.error_vel, delimiter)
        self.writeData(fname + ".dat", delimiter)
