import numpy
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.dates import MONDAY, SATURDAY
from pylab import meshgrid
import matplotlib.dates
import matplotlib.mathtext
import matplotlib.cm as cm
import matplotlib.axes as axes
from matplotlib.colors import LogNorm
from matplotlib import animation
import matplotlib.ticker 
import time, sys
from scipy.interpolate import interp1d
from utools.windrose import WindroseAxes
import utools.display_data as display_data
import utools.windows
from scipy.interpolate import interp1d
from matplotlib.patches import Rectangle
from scipy.odr.models import quadratic
import os, csv
# from numba import lowering

years = matplotlib.dates.YearLocator()  # every year
months = matplotlib.dates.MonthLocator()  # every month
yearsFmt = matplotlib.dates.DateFormatter('%Y')

# every monday
mondays = matplotlib.dates.WeekdayLocator(MONDAY)


debug = False
g_idx = 0

def draw_windrose(wd, ws, type, loc = 'best', fontsize = 10, unit = "[m/s]", angle = 67.5):
    def new_axes():
        fig = plt.figure(figsize = (8, 8), dpi = 80, facecolor = 'w', edgecolor = 'w')
        rect = [0.1, 0.1, 0.8, 0.8]
        wax = WindroseAxes(fig, rect, axisbg = 'w')
        fig.add_axes(wax)
        wax.set_rfontsize(fontsize-1)
        wax.rangle = angle
        return wax

    def set_legend(wax, loc, fontsize = 10, unit = "[m/s]"):
        l = wax.legend(loc = loc, unit = unit)  # axespad = -0.10)
        plt.setp(l.get_texts(), fontsize = fontsize)

    if type == 'bar':
        # windrose like a stacked histogram with normed (displayed in percent) results
        wax = new_axes()
        # for i in range(0, len(wd)):
        #    print "%d) direction: %f,  speed: %f " % (i, wd[i], ws[i])
        wax.bar(wd, ws, nsector = 32, normed = True, opening = 0.8, edgecolor = 'white')
        plt.xticks(fontsize = fontsize)
        set_legend(wax, loc, fontsize - 2, unit)

    if type == 'contour':
        # Same as above, but with contours over each filled region...
        wax = new_axes()
        wax.contourf(wd, ws, nsector = 32, bins = arange(0, 40, 4), cmap = cm.hot)
        wax.contour(wd, ws, nsector = 32, bins = arange(0, 40, 4), colors = 'black')
        set_legend(wax, loc, fontsize)

    if type == 'hist':
        wax = new_axes()
        wax.bar(wd, ws, normed = True, nsector = 32, opening = 0.8, edgecolor = 'white')
        set_legend(wax, loc, fontsize)

        table = wax._info['table']
        wd_freq = numpy.sum(table, axis = 0)
        dir = wax._info['dir']
        wd_freq = numpy.sum(table, axis = 0)

        # create another plain figure otherwise it will call the other bar() method from WindRose
        fig = plt.figure(figsize = (8, 8), dpi = 80, facecolor = 'w', edgecolor = 'w')
        ax = fig.add_subplot(111)
        ax.bar(numpy.arange(32), wd_freq, align = 'center')
        xlabels = ('N', '', 'NNE', '', 'N-E', '', 'ENE', '', 'E', '', 'ESE', '', 'S-E', '', 'SSE', '', 'S', '', 'SSW', '', 'S-W', '', 'WSW', '', 'W', '', 'WNW', '', 'N-W', '', 'NNW', '')
        xticks = numpy.arange(32)
        plt.gca().set_xticks(xticks)
        plt.draw()
        plt.gca().set_xticklabels(xlabels)
        plt.draw()

    # #print ax._info
    plt.show()
   
           
    return plt


def rotation_transform(u, v, tet, clockwise = False):
    '''
    Rotates counter clockwise (-tet) and counter clockwise (+tet)
    @param v: N-S velocity array + is toward N
    @param u: E-W velocity array + is toward E
    @param tet: angle in radians between N direction and the axis we want to project (e.g. longitudinal axis of a bay)

    @return: vp: the velocity along the axis of the channel \ ===> if counter clock wise (the oposit id clockwise)
             up: the velocity perpendicular on the axis     /   
    '''
    if clockwise:
        vp =  v * numpy.cos(tet) + u * numpy.sin(tet)    # ax Y
        up = -v * numpy.sin(tet) + u * numpy.cos(tet)   # ax X
        
        
    else:
        vp = v * numpy.cos(tet) - u * numpy.sin(tet)     # ax Y
        up = v * numpy.sin(tet) + u * numpy.cos(tet)   #as X
        
    return up, vp


def select_dates(fbegin, fend, dateTime, data):
    '''
        Select only records matching time from "begin" to "end"
        string dates are also converted
    '''
    sData = []
    sdateTime = []
    sz = len(data)
    for j in range(0, sz):
        srowData = []
        srowTime = []
        rdata = data[j]
        idx = 0
        for row in rdata:
            fdate = float(dateTime[idx])
            if fdate >= fbegin and fdate <= fend:
                try:
                    frow = float(row)
                    srowData.append(frow)
                    srowTime.append(fdate)
                except:
                    print("select_date error!")
            # endif
            idx += 1
        # end for
        sData.append(srowData)
        sdateTime.append(srowTime)
    # endfor
    return [sdateTime, sData]


def print_profile(name, profile, datetime, firstbin, Xrange= None, Yrange=None, interval=1, spline = False, save=False, dzdt = None):
    '''
    Prints an ADCP profile
    @param profile: array containing values at on time snapsot for all the bins
    @param fistbin: height form bottom of the first bin 
    @param interval: interval between bins 
    '''
    #shape: 'full', 'left', or 'right'
    global g_idx
    
    #shape ='full'
    #head_starts_at_zero=False
    #arrow_params={'length_includes_head':True, 'shape':shape, 'head_starts_at_zero':head_starts_at_zero}
    
    fig = plt.figure(facecolor = 'w', edgecolor = 'k')
    ax = fig.add_subplot(111)
    
    x = numpy.zeros(len(profile) + 1)
    y = numpy.zeros(len(profile) + 1)
    y[0] = 0
    
    #for i in range(0,len(profile)):
    for i in range(0,len(profile)):
        x[i+1] = profile[i]
        if i == 0:
            y[1] = firstbin
        else:  
            y[i+1] = y[i]+ interval
    
    if spline:
        itype = ["linear", "slinear",  "nearest", "cubic", "quadratic", "zero"]
        f = interp1d(y, x, kind=itype[1]) #slinear is the best option
        ynew = numpy.linspace(0, max(y), 40)
        ax.plot(f(ynew), ynew, linewidth = 0.3, color = 'r')
    else:
        ax.plot(x, y, linewidth = 0.3, color = 'r')

    ax.vlines([0], [0], [max(y)], linestyles='dashed', linewidths=2.0)
    #plot the arrows
    #print "****  Y size = %d , X size = %d" % (len(y), len(x))
    head_length = max(numpy.abs(x))/10.
    head_width = max(numpy.abs(y))/50.
    for i in range(0,len(profile)+1):
        if i != 0:
            #print "**idx=%d - %s -  Y[%d] = %f , X[%d]=%f **" % (g_idx, datetime, i, y[i], i, x[i])
            if abs(x[i]) <= 0.0000:
                x[i] = 0.0001
            ax.arrow( 0, y[i],  x[i],  0, fc="k", ec="k", lw =0.5, head_length=head_length, \
                      head_width=head_width, length_includes_head=True)
    
    ax.text(Xrange[0]+ 0.01, Yrange[1]-0.3 , datetime, fontsize=12)
    
    if dzdt != None:
        #plot the DZ in 1 h interval
        ax2 = ax.twinx()
        #ax2.bar([0], y[-1]-[dzdt], 0.04, color='b')
        ax2_ymax=0.6
        ypos = y[-1]*  ax2_ymax/Yrange[1]
        if dzdt > 0:
            facecolor ="blue"
        else:
            facecolor = "red"
            
        ax2.add_patch(Rectangle((0 - 0.02,  ypos ), 0.04, dzdt, facecolor=facecolor))
        ax2.set_ylim(0, ax2_ymax)
        ax2.set_ylabel("$\Delta$ Z [$m$]")    
 
    #Set Labels
    ax.set_xlabel("Velocity [$ms^{-1}$]")
    ax.set_ylabel("Depth [$m$]")
        
 
 
            
    if Yrange != None:
        ax.set_ylim(Yrange[0], Yrange[1])
    if Xrange != None:
        ax.set_xlim(Xrange[0], Xrange[1])
    if save:
        dattim = datetime.replace("/","-")
        dattim = dattim.replace(":","_")
        
        fname = name+"_" + dattim+".png"
        #fname = ("%s_%04d.png") %( name, g_idx)
        plt.savefig(fname)
        g_idx+=1
         #animateg git with convert
         #    convert -delay 20 -loop 0 *png animated.gif
         # or mp4
         #  ffmpeg -r 40 -i Cell3_%04d.png movie.mp4
    else:    
        plt.show()

def plot_temp_u_v(adcp, time1, evel, time2, nvel, time3, depth, interp = False):
    # 8) Mixed water, air ,img data
    custom = numpy.array(['Wind spd(m/s)', 'Air (T($^\circ$C)', 'Water T($^\circ$C)', ])
    # ToDO: Add short and long radiation
    print("Start display mixed subplots ")

    dateTimes1 = []
    data = []
    varnames = []
    ylabels1 = []
    dateTimes2 = []
    ylabels2 = []
    groups = []
    groupnames = []

    ylabels3 = ["Depth [m]", "Depth [m]", "Depth [m]"]
    revel = evel  # [::-1]
    rnvel = nvel  # [::-1]

    if adcp.goodbins == 0:
        lim1 = adcp.config.ranges[0]
        lim2 = adcp.config.ranges[len(adcp.config.ranges) - 1]
        lenght = len(adcp.config.ranges)
        imgs = [revel, rnvel, depth]
        dateTimes3 = [time1, time2 , time3]
    else:
        #=======================================================================
        # lim1 = adcp.config.ranges[0]
        # lim2 = adcp.config.ranges[adcp.goodbins] + .5
        # lenght = adcp.goodbins
        # imgs = [revel[:adcp.goodbins + 1], rnvel[:adcp.goodbins + 1] , depth]
        # dateTimes3 = [time1[:adcp.goodbins + 1], time2[:adcp.goodbins + 1], time3]
        #=======================================================================
        lim1 = adcp.config.ranges[0]
        lim2 = adcp.config.ranges[adcp.goodbins] - .5
        lenght = adcp.goodbins - 1
        imgs = [revel[:adcp.goodbins], rnvel[:adcp.goodbins] , depth]
        dateTimes3 = [time1[:adcp.goodbins], time2[:adcp.goodbins], time3]



    t11 = ['0', '2', '4', '6']
    t12 = [6, 4, 2, 0]
    t21 = ['0', '2', '4', '6']
    t22 = [6, 4, 2, 0]
    t31 = ['0', '3', '6', '9']
    t32 = [9, 6, 3, 0]
    tick = [[t11, t12], [t21, t22], [t31, t32]]



    maxdepth = [lim2, lim2, 9]
    mindepth = [lim1, lim1, 0]
    # revert the ADCP data


    firstlogdepth = [0, 0, 0]
    #===========================================================================
    # maxtemp = [0.3, 0.3 , 26]
    # mintemp = [-0.3, -0.3 , 0]
    #===========================================================================
    maxtemp = [0.2, 0.2 , 26]
    mintemp = [-0.2, -0.2 , 0]

    display_data.display_mixed_subplot(dateTimes1 = dateTimes1, data = data, varnames = varnames, ylabels1 = ylabels1,
                                       dateTimes2 = dateTimes2, groups = groups, groupnames = groupnames, ylabels2 = ylabels2,
                                       dateTimes3 = dateTimes3, imgs = imgs, ylabels3 = ylabels3, ticks = tick, maxdepths = maxdepth,
                                       firstlogs = firstlogdepth, maxtemps = maxtemp, mindepths = mindepth, mintemps = mintemp,
                                       fnames = None, interp = interp, revert = False, custom = None, maxdepth = None, tick = None, firstlog = None, yday = True, \
                                       title = False, grid = False, limits = None, sharex = True, fontsize = 20)

def plot_ADCP_velocity(adcp, data_args, graph_title, subplot = False, bins = None,
                       grid = False, title = False, sharex = True, doy = True) :
    '''UNTITLED Summary of this function goes here
      Detailed explanation goes here
    '''
    print('--> Once read this is what the ADCP data structure looks like!')

    # of graphs
    leng = int(len(data_args) / 2)

    # of data points
    L = numpy.size(adcp.number) - 1
    for i in range(0, leng):
        if subplot == 1:
            fig = plt.figure(num = i, facecolor = 'w', edgecolor = 'k')
        # end
        # plot(adcp.number, data_args(i).data);from pylab import *

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # plot against the velocity timeseries which are multiple dimension
        # matrix the first dimension is the data and the last the data
        # timeseries
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        data = data_args[0]
        label = data_args[1]
        if adcp.goodbins == 0:
            vel_bins = numpy.size(data[:, 0])
        else:
            vel_bins = numpy.size(data[:adcp.goodbins, 0])
            #vel_bins = numpy.size(data[:adcp.config.n_cells-1, 0])
        if (vel_bins > 9) and (subplot == 1):
            print('Too many graphs for subplot, changing to one graph per timeseries\n')
            subplot = 0
        # end
        timearr = adcp.mtime[0, 1:L - 1]
        if doy:
            dofy = numpy.zeros(len(timearr))
            for k in range(0, len(timearr)):
                d1 = matplotlib.dates.num2date(timearr[k])
                dofy[k] = d1.timetuple().tm_yday + d1.timetuple().tm_hour / 24. + d1.timetuple().tm_min / (24. * 60) + d1.timetuple().tm_sec / (24. * 3600)
            # end for
            x = dofy
        else :
            x = timearr
        # endif
        sp = 0
        for k in range(0, vel_bins):
            if bins is not None:
                if k not in bins:
                    continue
            if subplot == 1:
                if bins != None:
                    nob = len(bins)
                else:
                    nob = vel_bins
                subplt_num = 100 * (nob) + 10 + sp + 1
                ax = fig.add_subplot(subplt_num)
            else:
                fig = plt.figure(num = k, facecolor = 'w', edgecolor = 'r')
                ax = fig.add_subplot(111)
            # end
            y = data[k, 1:L - 1]
            ax.plot(x, y, linewidth = 0.3, color = 'r')

            if doy:
                ax.set_xlim(x[1], x[len(x) - 1])
            else:
                datemin = num2date(adcp.mtime[0, 1])
                datemax = num2date(adcp.mtime[0, L - 1])  # plotting range
                ax.set_xlim(datemin, datemax)

            if subplot == 1:
                ylabel = label + ' bin:%d' % k
            else:
                ylabel = label + ' bin:%d' % k
            plt.ylabel(ylabel)
            # plt.legend(ax, ylabel)

            if title and i == 0:
                plt.title(graph_title)
            # end

            if doy == False:
                # format the ticks
                formatter = matplotlib.dates.DateFormatter('%Y-%m-%d')
                # formatter = matplotlib.dates.DateFormatter('`%y')
                # ax.xaxis.set_major_locator(years)
                ax.xaxis.set_major_formatter(formatter)
                # ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%d'))
                ax.xaxis.set_minor_locator(mondays)

                # ax.xaxis.grid(True, 'major')
                ax.xaxis.grid(True, 'minor')
                # rotates and right aligns the x labels, and moves the bottom of the
                # axes up to make room for them
                fig.autofmt_xdate()

            if grid:
                ax.grid(True)


            # plt.show()
            sp += 1  # increment subplot.
        # end
        plt.show()

def  plot_ADCP_velocity_img(adcp, data_args, graph_title, echo = True, interp = False, doy = True
                            , sensor_height=0, water_depth=0):
        L = numpy.size(adcp.number) - 1
        #
        # PCOLOR
        #
        datemin = num2date(adcp.mtime[0, 1])

        done = False
        lim = L - 1
        while not done:
            try:
                datemax = num2date(adcp.mtime[0, lim])  # plotting range
                done = True
            except:
                lim = lim - 1

        data = numpy.mean(adcp.intens, 1)

        L = min(numpy.size(adcp.number[0]), numpy.size(data[0]))
        L = min(L, lim)

        if echo:
            displdata = data[:, :L - 1]  # this displays echo intensity
        else:
            displdata = data_args[0][:, :L - 1]  # This displays velocity

        displdata[numpy.isnan(displdata)] = 0  # set to zero NAN values
        displdata[numpy.isinf(displdata)] = 0  # set to 0 infinite values
        fig = plt.figure(100)
        ax = fig.add_subplot(111)

        timearr = adcp.mtime[0, 0:L - 1]

        if doy:
            dofy = numpy.zeros(len(timearr))
            for k in range(0, len(timearr)):
                d1 = matplotlib.dates.num2date(timearr[k])
                dofy[k] = d1.timetuple().tm_yday + d1.timetuple().tm_hour / 24. + d1.timetuple().tm_min / (24. * 60) + d1.timetuple().tm_sec / (24. * 3600)
            # end for
            x = dofy
        else :
            x = timearr

        # endif

        # select only the valid bins
        if adcp.goodbins == 0 or adcp.goodbins == -1:
            length = len(adcp.config.ranges)-1
            lim1 = adcp.config.ranges[0] + sensor_height
            lim2 = adcp.config.ranges[length] + sensor_height
            last = 0
            while lim2 > water_depth and water_depth != 0:
                length -= 1
                last -=1
                lim2 = adcp.config.ranges[length] + sensor_height
                displdata = data_args[0][:last, :L - 1]
        else:
            lim1 = adcp.config.ranges[0]
            lim2 = adcp.config.ranges[adcp.goodbins-1] + .5
            length = adcp.goodbins

        if interp == False:
            y = numpy.linspace(lim1, lim2, length)
            # X, Y = meshgrid(adcp.mtime[0][0:L - 1], adcp.config.ranges)
            X, Y = meshgrid(x, y)

            # displdata = data[0:L - 1][:]
            im = ax.pcolormesh(X, Y, displdata, shading = 'gouraud', vmin = -0.1, vmax = 0.1, cmap = cm.jet)
        else:
            multip = 5
            newdepth = numpy.linspace(lim1 - sensor_height, lim2 - sensor_height, length * multip)
            fint = interp1d(adcp.config.ranges[0:length+1], displdata.T, kind = 'cubic')
            newdata = fint(newdepth).T

            X, Y = meshgrid(x, newdepth)
            im = ax.pcolormesh(X, Y, newdata, shading = 'gouraud', vmin = -0.1, vmax = 0.1, cmap = cm.jet)

        cb = fig.colorbar(im)
        # cbar.ax.set_ylim([cbar.norm(3000), cbar.norm(6000)])
        # cbar.outline.set_ydata([cbar.norm(3000)] * 2 + [cbar.norm(6000)] * 4 + [cbar.norm(3000)] * 3)
        # cbar.ax.set_aspect(60)
        # cb.set_clim(-0.25, 0.25)
        # maxint = numpy.amax(displdata)
#===============================================================================
#         if echo:
#             cb.set_clim(0, maxint)
#         else:
#             cb.set_clim(-maxint / 4. , maxint / 4.)
#             # draw the v = 0
#             levels = [0.1]
#             colors = ['k']
#             linewidths = [0.5]
#             if interp == False:
#                 vel0 = displdata
#             else:
#                 vel0 = newdata
#
#             fontsize = 16
#             # ax.contour(X, Y, vel0, levels, colors = colors, linewidths = linewidths, fontsize = fontsize)
#===============================================================================


        labels = cb.ax.get_yticklabels()
        fontsize = 18
        plt.setp(labels, rotation = 0, fontsize = fontsize - 4)


        ylabel = ('Range (m)')

        plt.ylabel(ylabel).set_fontsize(fontsize + 2)
        title = ' Temperature Profiles'

        if debug == True:
            print("DATA")

            for j in range(0, 1000):
                print("%d" % j, end=' ')
                for i in range(0, data.shape[0]):
                    print(("i=%d : %f") % (i, data[i, j]), end=' ')
                print()

        if not doy:
            ax.xaxis_date()
            ax.set_xlim(datemin, datemax)

            # format the ticks
            formatter = matplotlib.dates.DateFormatter('%Y-%m-%d')
            ax.xaxis.set_major_formatter(formatter)
            ax.xaxis.set_minor_locator(mondays)
            ax.xaxis.grid(True, 'minor')
            ax.grid(True)

            # rotates and right aligns the x labels, and moves the bottom of the
            # axes up to make room for them
            fig.autofmt_xdate()


        plt.show()
    # end

def plotVelImg(fig, ax, axb, velarray,dateTime, maxdepth, interp=3, minval=None, maxval=None, revert =False,
               zerofirstline = True, colorbar = True, cblabel ="vel", fontsize =18):
    
    
    #insert zeros for thefirst roa in velarray 
    if zerofirstline:
        zeroline=numpy.zeros(len(velarray[0]))
        #first zero is the line number before to insert
        #The second zero is the axis
        #numpy.insert(velarray, 0, zeroline, 0) 
        velarray = numpy.vstack([zeroline, velarray])
    m = len(velarray)
    y = numpy.linspace(0, maxdepth, m)
    
    if interp:
        from scipy.interpolate import interp1d
        new_y = numpy.linspace(0, maxdepth, m * interp)
        fint = interp1d(y, velarray.T, kind = 'cubic')
        newVal = fint(new_y).T
        if revert == True:
            yrev = new_y[::-1]
        else:
            yrev = new_y
    else:
        if revert == True:
            yrev = y[::-1]
        else:
            yrev = y
        # if
        newVal = velarray
    # end if interp
    
    X, Y = numpy.meshgrid(dateTime, yrev)
    if maxval != None and minval != None:
        im = ax.pcolormesh(X, Y, newVal, shading = 'gouraud', vmin = minval, vmax = maxval)  # , cmap = 'gray', norm = LogNorm())
    else:
        im = ax.pcolormesh(X, Y, newVal, shading = 'gouraud')  # , cmap = 'gray', norm = LogNorm())

    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
    miny = numpy.min(Y)
    maxy = numpy.max(Y)
    dy = (maxy - miny)/4.0
    ax.set_yticks(numpy.arange(miny, maxy+dy/2., dy))

    if colorbar:
        cb = fig.colorbar(im, cax = axb)
        cb.set_clim(minval, maxval)
        labels = cb.ax.get_yticklabels()
        for t in labels:
            t.set_fontsize(fontsize - 8)
        from matplotlib import ticker
        tick_locator = ticker.MaxNLocator(nbins = 5)
        cb.locator = tick_locator
        cb.update_ticks()
        if cblabel:
            cb.set_label(cblabel)
            text = cb.ax.yaxis.label
            font = matplotlib.font_manager.FontProperties(size = fontsize -1)
            text.set_font_properties(font)
    # end if
   


def display_subplots(date, dateImg, dataarr, dnames = None, yday = None, tick = None, legend = None, hourgrid = False,
                     img=False, cbrange = [-0.8, 0.8], maxdepth = 5.0, fontsize = 16, minorlabel = False):
    '''
    @param date: one dimensional array containing the X - axis (time datapoints 
    '''

    matplotlib.mathtext.SHRINK_FACTOR = 0.9

    fig = plt.figure(facecolor = 'w', edgecolor = 'k')
    fformat = 100 * len(dataarr) + 10
    matplotlib.rcParams['legend.fancybox'] = True

    if yday :
        minx = 10000000.
        maxx = 0.

        # find maxX and minX of the X axis
        dmax = matplotlib.dates.num2date(date[len(date) - 1])
        maxx = max(maxx, (dmax.timetuple().tm_yday + dmax.timetuple().tm_hour / 24. + dmax.timetuple().tm_min / (24. * 60) + dmax.timetuple().tm_sec / (24. * 3600)))

        dmin = matplotlib.dates.num2date(date[0])
        minx = min(minx, (dmin.timetuple().tm_yday + dmin.timetuple().tm_hour / 24. + dmin.timetuple().tm_min / (24. * 60) + dmin.timetuple().tm_sec / (24. * 3600)))


    i = 0
    # ls = ['-', '--', ':', '-.', '-', '--', ':', '-.']
    ax = numpy.zeros(len(dataarr), dtype = axes.Subplot)

    for i in range(0,len(dataarr)):
        if i == 0:
            ax[i] = fig.add_subplot(fformat + i + 1)
        else:
            ax[i] = fig.add_subplot(fformat + i + 1, sharex=ax[0])
            
        depth = dataarr[i]
        # ax.plot(dateTime[1:], coef[1:])
        
        if yday:
            dofy = numpy.zeros(len(date))
            # dates = [datetime.fromordinal(d) for d in dataTime]
            # dofy = [d.tordinal() - datetime.date(d.year, 1, 1).toordinal() + 1 for d in dates]
            for j in range(0, len(date)) :
                d = matplotlib.dates.num2date(date[j])
                dofy[j] = d.timetuple().tm_yday + d.timetuple().tm_hour / 24. + d.timetuple().tm_min / (24. * 60) + \
                          d.timetuple().tm_sec / (24. * 3600)
            date =dofy
            
            dofyImg = numpy.zeros(len(dateImg))
            # dates = [datetime.fromordinal(d) for d in dataTime]
            # dofy = [d.tordinal() - datetime.date(d.year, 1, 1).toordinal() + 1 for d in dates]
            for j in range(0, len(dateImg)) :
                d = matplotlib.dates.num2date(dateImg[j])
                dofyImg[j] = d.timetuple().tm_yday + d.timetuple().tm_hour / 24. + d.timetuple().tm_min / (24. * 60) + \
                             d.timetuple().tm_sec / (24. * 3600)
            dateImg =dofyImg

        #-------------------------------------------------------------
        # X-AXIS -Time
        # format the ticks
        if yday == False:
            if hourgrid:
                formatter = matplotlib.dates.DateFormatter('%H')
            else:
                formatter = matplotlib.dates.DateFormatter('%Y-%m-%d')
            
            # formatter = dates.DateFormatter('`%y')
            ax[i].xaxis.set_major_formatter(formatter)
            # ax.xaxis.set_minor_formatter(dates.DateFormatter('%d'))
            #ax[i].xaxis.set_minor_locator(mondays)
            ax[i].xaxis.set_major_locator(matplotlib.dates.HourLocator())
        else:
            if hourgrid:
                majorLocator   = matplotlib.dates.DayLocator()
                majorFormatter = matplotlib.ticker.FormatStrFormatter('%d')
                minorLocator   = matplotlib.dates.HourLocator() 
                minorFormatter = matplotlib.dates.DateFormatter('%H')
                ax[i].xaxis.set_minor_locator(minorLocator)
                if minorlabel:
                    ax[i].xaxis.set_minor_formatter(minorFormatter)
            else:
                majorLocator   = matplotlib.dates.DayLocator() 
                majorFormatter = matplotlib.ticker.FormatStrFormatter('%.02f')
                minorLocator   = matplotlib.dates.DayLocator() 
                minorFormatter = matplotlib.ticker.FormatStrFormatter('%.02f')
            
            ax[i].xaxis.set_major_locator(majorLocator)
            ax[i].xaxis.set_major_formatter(majorFormatter)
            
        if dnames != None:
            ylabel = dnames[i]
        else:
            ylabel = 'Depth. (m)'
            
        ax[i].set_ylabel(ylabel).set_fontsize(fontsize)
       
        if yday == False:
            minx = date[0]
            maxx = date[len(date) - 1]
            ax[i].set_xlim(xmin=minx, xmax = maxx)
            if i == len(ax)-1:
                ax[i].set_xlabel("Time [hours]").set_fontsize(fontsize)
            #ax[i].set_xlabel("Time").set_fontsize(20)
            dx = (maxx+2 - minx)/23.0
            lower = minx
            upper =  maxx+dx/50.
            
            ax[i].set_xlim(xmin = lower, xmax = upper)
            ax[i].set_xticks(numpy.arange(lower, upper, dx))
            
        else:
            if i == len(ax)-1:
                ax[i].set_xlabel("day of the year").set_fontsize(fontsize)

        
        # ax[i].legend(lplt, title = lg, shadow = True, fancybox = True)
        if tick != None:
            ax[i].set_yticks(tick[i][1])
            ax[i].set_yticklabels(tick[i][0])
       
        ax[i].xaxis.grid(which='minor', alpha=0.3)                                                
        ax[i].xaxis.grid(which='major', alpha=0.6)    
        
        if i != len(dataarr)-1:
            plt.setp( ax[i].get_xticklabels(), visible=False)
        #ax[i].xaxis.set_minor_locator(matplotlib.dates.HourLocator() )
        
        #--------------------------------------------------------------
        for tk in ax[i].xaxis.get_major_ticks():
                tk.label.set_fontsize(fontsize - 1)
        for tk in ax[i].yaxis.get_major_ticks():
                tk.label.set_fontsize(fontsize - 1)
        

        if type(depth) is list:
            ndepth = numpy.array(depth)
        else:
            ndepth = depth  
        if ndepth.ndim > 1:
            if img:
                # Make some room for the colorbar
                fig.subplots_adjust(left = 0.1, right = 0.86)

                # Add the colorbar outside...
                box = ax[i].get_position()
                pad, width = 0.016, 0.015
                axb = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
                plotVelImg(fig, ax[i], axb, ndepth, dateImg, maxdepth=maxdepth, interp=4, minval=cbrange[0],
                           maxval=cbrange[1], zerofirstline = True, colorbar = True, cblabel=r"$\mathsf{V_{alg} [m\ s^{-1}}$]",
                           fontsize = fontsize)
            else:
                for j in range(0, len(ndepth)):
                    ax[i].plot(dateImg, ndepth[j], color='k', linewidth = 1.4)
                    miny = numpy.min(ndepth[j])
                    maxy = numpy.max(ndepth[j])
                    ax[i].yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.3f'))
                    dy = (maxy - miny)/4.0
                    ax[i].set_yticks(numpy.arange(miny, maxy+dy/2., dy))
                    
        else:
            ax[i].plot(date, ndepth, color='k', linewidth = 1.4)
            if i == 0:
                ax[i].yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.2f'))
            else:
                ax[i].yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.3f'))
            miny = numpy.min(ndepth)
            maxy = numpy.max(ndepth)
            dy = (maxy - miny)/4.0
            ax[i].set_yticks(numpy.arange(miny, maxy+dy/2., dy))

        # LEGEND
        if not img:
            if type(legend[i]) is list:
                text = legend[i]
            else:
                text = [legend[i]] 
            ax[i].legend(text, loc = 'best', shadow = True, fancybox = True)


#===============================================================================
#         # X-AXIS -Time
#         # format the ticks
#         if yday == False:
#             if hourgrid:
#                 formatter = matplotlib.dates.DateFormatter('%H')
#             else:
#                 formatter = matplotlib.dates.DateFormatter('%Y-%m-%d')
#             
#             # formatter = dates.DateFormatter('`%y')
# 
#             ax[i].xaxis.set_major_formatter(formatter)
#             # ax.xaxis.set_minor_formatter(dates.DateFormatter('%d'))
#             #ax[i].xaxis.set_minor_locator(mondays)
#             ax[i].xaxis.set_major_locator(matplotlib.dates.HourLocator())
#         else:
#             if hourgrid:
#                 majorLocator   = matplotlib.dates.DayLocator()
#                 majorFormatter = matplotlib.ticker.FormatStrFormatter('%d')
#                 minorLocator   = matplotlib.dates.HourLocator() 
#                 minorFormatter = matplotlib.dates.DateFormatter('%H')
#                 ax[i].xaxis.set_minor_locator(minorLocator)
#                 if minorlabel:
#                     ax[i].xaxis.set_minor_formatter(minorFormatter)
#             else:
#                 majorLocator   = matplotlib.dates.DayLocator() 
#                 majorFormatter = matplotlib.ticker.FormatStrFormatter('%.02f')
#                 minorLocator   = matplotlib.dates.DayLocator() 
#                 minorFormatter = matplotlib.ticker.FormatStrFormatter('%.02f')
#             
#             ax[i].xaxis.set_major_locator(majorLocator)
#             ax[i].xaxis.set_major_formatter(majorFormatter)
#             
#         if dnames != None:
#             ylabel = dnames[i]
#         else:
#             ylabel = 'Depth. (m)'
#             
#         ax[i].set_ylabel(ylabel).set_fontsize(fontsize)
#        
#         if yday == False:
#             minx = date[0]
#             maxx = date[len(date) - 1]
#             ax[i].set_xlim(xmin=minx, xmax = maxx)
#             if i == len(ax)-1:
#                 ax[i].set_xlabel("Time [$hours$]").set_fontsize(fontsize)
#             #ax[i].set_xlabel("Time").set_fontsize(20)
#             ax[i].set_xlim(xmin = minx, xmax = maxx)
#             dx = (maxx - minx)/6.0
#             ax[i].set_xticks(numpy.arange(minx, maxx+dx/2., dx))
#             
#         else:
#             if i == len(ax)-1:
#                 ax[i].set_xlabel("day of the year").set_fontsize(fontsize)
# 
#         
#         # ax[i].legend(lplt, title = lg, shadow = True, fancybox = True)
#         if tick != None:
#             ax[i].set_yticks(tick[i][1])
#             ax[i].set_yticklabels(tick[i][0])
#        
#        
#         ax[i].xaxis.grid(which='minor', alpha=0.3)                                                
#         ax[i].xaxis.grid(which='major', alpha=0.6)    
#         
#         if i != len(dataarr)-1:
#             plt.setp( ax[i].get_xticklabels(), visible=False)
#         #ax[i].xaxis.set_minor_locator(matplotlib.dates.HourLocator() )
#         
#         plt.setp(ax[i].get_xticklabels(), visible = True, fontsize = fontsize-3)
#         plt.setp(ax[i].get_yticklabels(), visible = True, fontsize = fontsize-3)
#===============================================================================
        i += 1

    # end for

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them

    #if yday == None:
    #    fig.autofmt_xdate()
    plt.xticks(fontsize = fontsize - 1)
    plt.yticks(fontsize = fontsize - 1)
    plt.draw()
    plt.show()

def writeVeltoCSV(fname, data, append = False):
    if append:
        ofile = open(fname, "at")
    else:
        ofile = open(fname, "wt")
    writer = csv.writer(ofile, delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
    writer.writerow(data)
    ofile.close()

def ReadVolfromCSV(fname, level):
    ifile = open(fname, "rt")
    reader = csv.reader(ifile, delimiter = ',', quotechar = '"')
    for row in reader:
        if int(row[0]) == level:
            maxdepth = float(row[1])
            volume = float(row[2])
            area = float(row[3])
            break
    ifile.close()
    return maxdepth, volume, area
