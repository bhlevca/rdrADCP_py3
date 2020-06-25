import numpy
import matplotlib.pyplot as plt
import matplotlib.axes as axes
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.dates import MONDAY, SATURDAY
from pylab import meshgrid
import matplotlib.dates
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

import time
import math, sys
from scipy import fftpack



import wavelets.kCwt
import ufft.FFTGraphs as fftGraphs

from . import readRawBinADCP
debug = False


def nextpow2(i):
    """
    Find the next power of two

    >>> nextpow2(5)
    8
    >>> nextpow2(250)
    256
    """
    # do not use numpy here, math is much faster for single values
    buf = math.ceil(math.log(i) / math.log(2))
    return int(math.pow(2, buf))

'''
def nextpow2(i):
    """
    Find 2^n that is equal to or greater than.
    """
    n = 2
    while n < i: n = n * 2
    return n
'''

def plot_velocity_wavelet_spectrum(time, wdata, scaleunit = 'day'):
    kwavelet = wavelets.kCwt.kCwt(time, wdata, scaleunit, False)

    dj = 0.05  # Four sub-octaves per octaves
    s0 = -1  # 2 * dt                              # Starting scale, here 6 months
    J = -1  # 7 / dj                               # Seven powers of two with dj sub-octaves
    alpha = 0.5  # Lag-1 autocorrelation for white noise
    slevel = 0.95
    val1, val2 = (0, 260000)  # Range of sc_type (ex periods) to plot in spectogram
    avg1, avg2 = (0, 260000)  # Range of sc_type (ex periods) to plot in spectogram

    # [wave, scales, freq, coi, fft, fftfreqs, iwave, power, fft_power, amplitude, phase] = \
    kwavelet.doSpectralAnalysis('Velocity Wavelet Spectrum', "morlet", slevel, avg1, avg2, dj, s0, J, alpha, False)

    ylabel_ts = "amplitude"
    yunits_ts = 'm/s'
    xlabel_sc = ""
    ylabel_sc = 'Period [%s]' % kwavelet.wpar1.tunits
    # ylabel_sc = 'Freq (Hz)'
    sc_type = "period"
    # sc_type = "freq"
    # x_type = 'date'
    x_type = 'dayofyear'

    # we want to have the vector, power is in absolute value
    kwavelet.plotSpectrogram(kwavelet.wpar1, ylabel_ts, yunits_ts, xlabel_sc, ylabel_sc, sc_type, x_type, val1, val2, \
                             raw = False, title = False, powerimg = False, bcolbar = False)

    # ylabel_sc = 'Frequency ($s^{-1}$)'
    # kwavelet.plotAmplitudeSpectrogram(ylabel_ts, yunits_ts, xlabel_sc, ylabel_sc, sc_type, x_type, val1, val2)



def plot_rotary_wavelet_spectrum(adcp, bin = bin, counterclockwise = True, scaleunit = 'day'):
    print('--> Depth averaged velocity analysis')
    L = numpy.size(adcp.number[0]) - 1 - int(numpy.size(adcp.number[0]) / 10)

    t = adcp.mtime[0][0:L]
    u = numpy.nan_to_num(adcp.east_vel[bin])[0:L]
    v = numpy.nan_to_num(adcp.north_vel[bin])[0:L]

    kwavelet = wavelets.kCwt.kCwt(t, u, v, scaleunit, False)
    dj = 0.05  # Four sub-octaves per octaves
    s0 = -1  # 2 * dt                              # Starting scale, here 6 months
    J = -1  # 7 / dj                               # Seven powers of two with dj sub-octaves
    alpha = 0.5  # Lag-1 autocorrelation for white noise
    slevel = 0.95
    val1, val2 = (0, 260000)  # Range of sc_type (ex periods) to plot in spectogram
    avg1, avg2 = (0, 260000)  # Range of sc_type (ex periods) to plot in spectogram
    if counterclockwise:
        title = "Counterclockwise velocity spectrum"
    else:
        title = "Clockwise velocity spectrum"

    [wave, scales, freq, coi, fft, fftfreqs, iwave, power, fft_power, amplitude, phase] = \
        kwavelet.doSpectralAnalysis(title, "morlet", slevel, avg1, avg2, dj, s0, J, alpha, counterclockwise)

    ylabel_ts = "amplitude"
    yunits_ts = 'm/s'
    xlabel_sc = ""
    ylabel_sc = 'Period [%s]' % kwavelet.wpar1.tunits
    # ylabel_sc = 'Freq (Hz)'
    sc_type = "period"
    # sc_type = "freq"
    # x_type = 'date'
    x_type = 'dayofyear'

    # we want to have the vector, power is in absolute value
    kwavelet.plotSpectrogram(ylabel_ts, yunits_ts, xlabel_sc, ylabel_sc, sc_type, x_type, val1, val2, powerimg = False)
    # ylabel_sc = 'Frequency ($s^{-1}$)'
    # kwavelet.plotAmplitudeSpectrogram(ylabel_ts, yunits_ts, xlabel_sc, ylabel_sc, sc_type, x_type, val1, val2)


def plot_cross_spectogram_w_T(t, u, v, tt, T, scaleunit = 'day'):
    '''
    Plot the  cross wavelet spectrogram between the ADCP velocity and Temperature
    :param t : time array
    :param u : E velocity
    :param v : N velocity
    :param T :
    '''

    w = u + 1j * v

    if len(w) != len(T):
        ratio = len(w) / len(T)
        print("len w: %d | len t: %d | len tt: %d | len T: %d" % (len(w), len (t), len(tt), len(T)))
        print("ratio: %f" % ratio)
        oldrange = list(range(0, ratio * len(t), ratio))
        print("len(oldrange) :%d" % len(oldrange))
        iT = numpy.interp(t, tt, T)
        # iT2 = numpy.interp(range(0, ratio * len(t)), oldrange), T)
        print("len w: %d | len iT: %d |  len T: %d" % (len(w), len(iT), len(T)))
        tt = t
    else:
        iT = T

    kwavelet = wavelets.kCwt.kCwt(t, w, tt, iT, scaleunit, False, True)
    dj = 0.05  # Four sub-octaves per octaves
    s0 = -1  # 2 * dt                              # Starting scale, here 6 months
    J = -1  # 7 / dj                               # Seven powers of two with dj sub-octaves
    alpha = 0.5  # Lag-1 autocorrelation for white noise
    slevel = 0.95
    val1, val2 = (0, 520000)  # Range of sc_type (ex periods) to plot in spectogram
    avg1, avg2 = (0, 520000)  # Range of sc_type (ex periods) to plot in spectogram
    title = "Cross vel-T spectrum"
    kwavelet.doCrossSpectralAnalysis(title, "morlet", slevel, avg1, avg2, dj, s0, J, alpha)
    ylabel_ts = "PSD"
    yunits_ts = 'm/s'
    xlabel_sc = ""
    ylabel_sc = 'Period [%s]' % kwavelet.wpar1.tunits
    sc_type = "period"
    # sc_type = "freq"
    # x_type = 'date'
    x_type = 'dayofyear'


    print("#plot cross wavelet spectrogram")
    kwavelet.plotXSpectrogram(kwavelet.get_xwt(), extend = 'both', da = [6, 600])

    print("plot coherence wavelet spectrogram")
    kwavelet.plotXSpectrogram(kwavelet.get_wct(), extend = 'neither', x_type = x_type, ylabel_sc = ylabel_sc, da = [6, 600],
                  tfactor = kwavelet.wpar1.tfactor, crange = numpy.arange(0, 1.1, 0.1), scale = 'linear', angle = kwavelet.get_wct().angle)



def plot_depth_averaged_analysis(adcp, data_args, graph_title, subplot = False, scale = 'log', bin = 6, avg = False) :
    '''UNTITLED Summary of this function goes here
      Detailed explanation goes here
    '''
    print('--> Depth averaged velocity analysis')
    vel, ylabel = data_args
    L = numpy.size(adcp.number[0]) - 1

    if avg:
        averaged_vel = readRawBinADCP.nmean(vel, 0)
    else:
        averaged_vel = numpy.nan_to_num(vel[bin])


    t = adcp.mtime[0][0:L]
    averaged_vel = averaged_vel[0:L]

    # show extended calculation of spectrum analysis
    plot_fft_analysis(t, averaged_vel, graph_title, ylabel, subplot = subplot, scale = scale, bin = bin, avg = avg)


def plot_fft_analysis(t, ydata, graph_title, ylabel, subplot = False, scale = 'log', bin = 6, avg = False) :
    show = True
    data = [t, ydata]
    fftsa = fftGraphs.FFTGraphs(None, None, None, show, data)
    showLevels = False
    detrend = False
    draw = False
    tunits = "day"
    funits = "cph"
    filter = None
    window = "hanning"
    num_segments = 4

    [Time, y, x05, x95, fftx, freq, mx] = fftsa.doSpectralAnalysis(showLevels, draw, tunits, window, num_segments, filter, scale)
    fftsa.plotSingleSideAplitudeSpectrumFreq(ylabel, None, funits, y_label = ylabel, title = graph_title, log = scale, \
                                            fontsize = 20, tunits = tunits, plottitle = True, grid = True, \
                                            ymax = None, graph = True)

#===============================================================================
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     # ax.plot(tph, mx)
#     # plt.xlabel("Period (h)")
#     # ax.plot(f * 3600, mx)
#     ax.set_yscale(scale)
#
#     ax.loglog(f * 3600, mx, basex = 10)
#     plt.xlabel("Frequency [cph]")
#     plt.ylabel(ylabel)
#     ax.xaxis. set_major_locator(MaxNLocator(12))
#     ax.xaxis.grid(True, 'minor')
#     ax.grid(True)
#     plt.show()
#===============================================================================
# end


def plot_FFT_twinx_W_T(time, W, Ttime, T, scale = 'log', drawslope = False) :
    '''UNTITLED Summary of this function goes here
      Detailed explanation goes here
    '''
    show = True
    draw = False

    if len(W) != len(T):
        ratio = len(W) / len(T)
        print("len w: %d | len t: %d | len tt: %d | len T: %d" % (len(W), len (time), len(Ttime), len(T)))
        print("ratio: %f" % ratio)
        iT = numpy.interp(time, Ttime, T)
        print("len w: %d | len iT: %d |  len T: %d" % (len(W), len(iT), len(T)))
    else:
        iT = T

    data = [time, W]
    data1 = [time, iT]

    fftsa = fftGraphs.FFTGraphs(None, None, None, show, data, data1)
    showLevels = False
    tunits = "day"
    funits = "cph"
    log = True
    filter = None
    window = "hanning"
    num_segments = 10

    [Time, y, x05, x95, fftx, freq, mx] = fftsa.doSpectralAnalysis(showLevels, draw, tunits, window, num_segments, filter, log)
    fftsa.plotPSD2SpectrumsFreq("Velocity", "Temperature", funits = funits, title = None,
                                         y_label1 = "Velocity PSD [m$^2$/s]", tunits = 'day',
                                         log = log, fontsize = 20, plottitle = False, grid = True, \
                                         ymax = None, graph = True,
                                         twoaxes = True, ylabel2 = 'Temperature PSD[$\circ$C$^2$/s]', ymax_lim2 = None, drawslope = drawslope)
# end


def plot_FFT_Three_V_T_WL(time, V, Ttime, T, WL, WTime, scale = 'log', drawslope = False, plotTemp=True) :
    '''UNTITLED Summary of this function goes here
      Detailed explanation goes here
    '''
    show = True
    draw = False

    if len(V) != len(T):
        ratio = len(V) / len(T)
        print("len w: %d | len t: %d | len tt: %d | len T: %d" % (len(V), len (time), len(Ttime), len(T)))
        print("ratio: %f" % ratio)
        iT = numpy.interp(time, Ttime, T)
        print("len w: %d | len iT: %d |  len T: %d" % (len(V), len(iT), len(T)))
    else:
        iT = T

    if len(V) != len(WL):
        ratio = len(V) / len(WL)
        print("len w: %d | len t: %d | len tt: %d | len T: %d" % (len(V), len (time), len(WTime), len(WL)))
        print("ratio: %f" % ratio)
        iWL = numpy.interp(time, WTime, WL)
        print("len w: %d | len iT: %d |  len T: %d" % (len(V), len(iT), len(T)))
    else:
        iWL = WL


    data = [time, V]
    data1 = [time, iWL]
    data2 = [time, iT]

    fftsa = fftGraphs.FFTGraphs(None, None, None, show, data, data1, data2)
    showLevels = False
    tunits = "day"
    funits = "cph"
    filter = None
    window = "hanning"
    num_segments = 10
    grid = False

    [Time, y, x05, x95, fftx, freq, mx] = fftsa.doSpectralAnalysis(showLevels, draw, tunits, window, num_segments, filter, scale)
    #===========================================================================
    # This is for PSD
    # fftsa.plotPSD3SpectrumsFreq("Velocity", "Temperature", "Water Level", funits = funits, \
    #                                      log = scale, fontsize = 20, plottitle = False, grid = True, \
    #                                      ymax = None, graph = True, tunits = tunits, ymax_lim2 = None, \
    #                                      ylabel1 = "Velocity PSD [m$^4 s^{-3}$]", \
    #                                      ylabel2 = 'Temperature PSD[$\circ$C$^2 s^{-1}$]', \
    #                                      ylabel3 = "Water level PSD [m$^2 s^{-1}$]", 
    #                                      ymax_lim3 = None, drawslope = drawslope)
    #===========================================================================
    
    # this is for single side amplitude FFT
    if plotTemp:
        fftsa.plotPSD3SpectrumsFreq("Velocity", "Temperature", "Water Level", funits = funits, \
                                         log = scale, fontsize = 20, plottitle = False, grid = grid, \
                                         ymax = None, graph = True, tunits=tunits, ymax_lim2 = None, \
                                         ylabel1 = "Velocity [ms$^{-1}$]", \
                                         ylabel2 = "Water level [m]", \
                                         ylabel3 = 'Temperature [$^{\circ}$C]', \
                                         ymax_lim3 = None, drawslope = drawslope)
    else:
        fftsa.plotPSD2SpectrumsFreq("Velocity", "Water Level", funits=funits,
                                         ylabel1="Velocity [ms$^{-1}$]",
                                         ylabel2="Water level [m]",
                                         log=scale, fontsize = 20, plottitle=False, grid=grid, \
                                         ymax=None, graph=True,
                                         ymax_lim2=None, drawslope=drawslope)
# end


